function [d, S] = solve(sol, param, solver_type, tem_info, trafo)
    % Solve TEM fwp and provides sensitivity matrix.
    %
    % Sensitivity matrix = derivative of observed quantity (i.e. apparent
    %                      resistivities) w.r.t. cell parameters.
    %
    % SYNTAX
    %   [d, S] = solve(sol, type, dc_info, param, trafo)
    %
    % INPUT PARAMETER
    %   sol ... Struct, containing
    %              O          - Matrix (observation operator).
    %              dofmap     - Struct, containing cell2DOF mapping as well
    %                           as the mesh and FE element objects.
    %              type       - Character, denoting problem dimension / type.
    %              KMB_handle - Handle taking vector 'param' providing
    %                           system matrix A and rhs vector b
    %              TM         - Mass matrix derivative tensor
    %   param   ... Vector of cell parameters.
    %HA
    % OPTIONAL PARAMETER
    %   solver_type ... Character, denoting the solver for linear system.
    %   tem_info    ... Struct, containing survey info and mappings to obtain
    %                   observation from simulated quantities.
    %   trafo       ... Handle, taking matrix and vector param and applies
    %                   transformation of param on matrix.
    %
    % OUTPUT PARAMETER
    %   d ... Solution vector [n_data, 1] with n_data = n_dim*n_RX
    %   S ... Sensitivity matrix [n_data, n_cells].
    %
    % REMARKS
    %   sol      is provided by app_tem.fwd.assemble.m
    %   tem_info is provided by app_tem.data.get_obs.m
    %
    % The elements of the data vector d are sorted as follows:
    %    vector (field) components, i.e. dofmap.mesh.dim    ... i
    %    observation points                                 ... k
    %    data types                                         ... l
    %    observation time steps                             ... t
    %
    %    d_i,k,l,t
    %    i.e. reshape(d, n_dim, n_RX, n_data_type, n_t) provides sorted
    %    4-way tensor of data where each quantity can easily be accessed.

    %% Fetch info / set defaults and check input.

    assert(isstruct(sol) && all(isfield(sol, {'KMb_handle', 'type'})) ...
           && any(strcmp(sol.type, {'3D'})));
    if ~isfield(sol, 'O') && nargout == 2
        error('Field sol.O required to provide sensitivity matrix.');
    end
    assert(isstruct(tem_info) && all(isfield(tem_info, {'t_obs', 't', 'pt_TX'})));

    if isempty(solver_type)
    % -> linear system solver
       solver_type = 'backslash';
    else
       assert(ischar(solver_type));
    end

    if ~isfield(tem_info, 'integration_type') || isempty(tem_info.integration_type)
    % -> time integration approach
        tem_info.integration_type = 'impl';
    else
        assert(ischar(tem_info.integration_type));
    end

    % Sanity check.
    all(ismember(tem_info.t_obs, tem_info.t));

    % Set up parallel pool.
    % Note: only for sensitivity calculation parallelism can be exploited
    if nargout > 1
        if ~isfield(tem_info, 'solver_parallel') || tem_info.solver_parallel <= 1
            % No parallelism required.
            tem_info.solver_parallel = 0;
        elseif tem_info.solver_parallel > 1
            % Check if parallel toolbox license is available.
            if license('test', 'distrib_computing_toolbox')
                % Check for existing pool of parallel worker.
                par_pool = gcp('nocreate');
                if isempty(par_pool)
                    % Create pool.
                    par_pool = parpool('local', tem_info.solver_parallel);
                elseif (par_pool.NumWorkers == tem_info.solver_parallel)
                    fprintf(['Using active pool with %d worker. ', ...
                            'Use delete(gcp(''nocreate'')); to close current pool.\n'], par_pool.NumWorkers);
                else
                    delete(par_pool);
                    par_pool = parpool('local', tem_info.solver_parallel);
                end
                % Ensure that during (back)propagation via parfor only
                % single indices (and no chunks) are passed, such that
                % workload between the worker can be better balnaced
                parpool_opts = parforOptions(par_pool, 'RangePartitionMethod', 'fixed', 'SubrangeSize', 1);
            else
                warning(['No Distributed Computing Toolbox available. ', ...
                        'Proceed with serial calculation.']);
                tem_info.solver_parallel = 0;
            end
        end
    end

    %% Get linear system info.

    % Fetch system and measurement operator.
    [K, M, b] = sol.KMb_handle(param);
    O = sol.O;

    %% Handle time integration.

    switch tem_info.integration_type
        case 'impl'
            % Implicit euler.
            if ~isfield(tem_info, 'src_info')
                src_info = struct();
                src_info.type = 'shut-off';         % default rhs function
                src_info.dt_ref = tem_info.t(1)/10; % default time refinement close to t=0
                tem_info.src_info = src_info;
            else
                assert(isstruct(tem_info.src_info));
                assert(all(isfield(tem_info.src_info, {'type', 'dt_ref'})));
            end
            % Check if corresponding dt is given.
            % Note: Calculating dt from t (e.g. using diff()) is not
            %       recommended as due to rounding errors nearly each time
            %       step will be associated with a unique dt.
            assert(isfield(tem_info, 'dt'));
            assert(isvector(tem_info.dt) && ...
                   length(tem_info.dt) == length(tem_info.t));

            % Solve fwp.
            if nargout > 1
                [d, S] = solve_impl(b, K, M, O, sol.TM, solver_type, tem_info);

                % Apply chain rule for derivative of parameter w.r.t. a trafo.
                if nargin == 5
                    S = trafo(S, param);
                end
            else
                d = solve_impl(b, K, M, O, sol.TM, solver_type, tem_info);
            end

        case 'rba'
            % Rational best approximation.
            % Note: only shut-off rhs can be represented.
            if isfield(tem_info, 'n_poles')
                assert(isscalar(tem_info.n_poles));
            else
                % Set default number of poles
                % Note: Too less poles may result in an artificial field at
                %       the source position at late times!
                tem_info.n_poles = 12;
            end

            % Solve fwp.
            if nargout > 1
                [d, S] = solve_rba(b, K, M, O, solver_type, tem_info);

                % Apply chain rule for derivative of parameter w.r.t. a trafo.
                if nargin == 5
                    S = trafo(S, param);
                end
            else
                d = solve_rba(b, K, M, O, solver_type, tem_info);
            end

        otherwise
            error('Unknown time integration approach.');
    end

    %% Helper.

    function solver = get_solver_fun(A, type, num)
        % Wrapper for solver classes.

        % Define uniform solver handle interface.
        switch type
            case {'backslash', 'chol'}
                if strcmp(type, 'backslash')
                    type_ = 'auto';
                else
                    type_ = type;
                end
                solver_ = decomposition(A, type_);
                solver = struct();
                solver.solve = @(b) solver_ \ b;
                solver.delete = @() []; % dummy, .solve will not be touched!

            case  {'mumps'}
                solver = solving.MUMPS(A);

            case  {'mumps_ooc'}
                solver = solving.MUMPS_OOC(A, num);

            otherwise
                error('Unsupported solver type.');
        end
    end

    function [d, S] = solve_rba(b, K, M, O, solver_type, tem_info)
	error('Not implemented, yet!')
    end

    function [d, S] = solve_impl(b, K, M, O, TM, solver_type, tem_info)
        % Create forward solution and sensitivity matrix.
        %
        % t = 0 is assumed as initial time

        % Get time-dependent rhs.
        b_fun = get_djdt_fun(tem_info.src_info);

        % Refine temporal resolution around t = 0:
        [t, dt] = adjust_t_dt(tem_info);

        % Fetch infos.
        dt_unique = unique(dt);
        n_dt = numel(dt_unique);
        n_t = numel(t);
        n_t_obs = numel(tem_info.t_obs);
        n_param = size(TM, 3);
        n_dof = size(K, 1);
        n_d = size(O, 2);

        % Get indices of fwd calculation time steps that correspond to t_obs.
        [tmp, t_obs_idx_in_t] = ismember(t, tem_info.t_obs);
        assert(all(sort(tem_info.t_obs) == tem_info.t_obs)); % assume increasing times
        assert(numel(find(tmp)) == length(tem_info.t_obs));

        % Assign unique solver, i.e. system decomposotions, to all time steps.
        solver_idx_in_t = dt.' == dt_unique;
        solver_idx_in_t = (solver_idx_in_t * (1:n_dt).').';
        assert(length(solver_idx_in_t) == n_t, ...
               'Ensure dt to have same size as t.');

        % Create all unique solver objects (if ooc, additionally store
        % factorizations to disc).
        solver_MdtK = cell(n_dt, 1);
        fprintf('Calculate %d system decompositions for %d DOF: ', n_dt, size(M, 1));
        run_time = tic;
        for dd = 1:n_dt
            str_verbose = fprintf('%d', dd);
            % Note: Variable dd only relevant for ooc (denotes file name),
            %       otherwise it's just a dummy argument.
            solver_MdtK{dd} = get_solver_fun(dt_unique(dd)*K + M, solver_type, dd);
            fprintf(repmat('\b', 1, str_verbose));
        end
        fprintf('%.2ds \n', toc(run_time));

        %%% Calculate data via forward propagation %%%

        % Loop over all time steps.
        u = zeros(n_dof, n_t_obs);
        d_ = zeros(n_d, n_t_obs);
        fprintf('Propagate primal solution %d time steps w.r.t. %d RHS: ', n_t, size(b, 2));
        run_time = tic;
        for tt = 1:n_t

            % Solve linear system arising from implicit Euler scheme.
            % ->   u_1 = (M + dt*K)^-1 * [(M*u0) + dt*dj/dt]
            %    M*u_0 = j(t=0) + M*u_DC
            str_verbose = fprintf('%d', tt);
            if tt == 1
                % Initial time step:
                if isfield(tem_info, 'u_DC')
                    % Incorporate an electrical field distribution that
                    % exists for t < 0
                    % -> M*u_0 = ... + M*u_DC
                    assert(isvector(tem_info.u_DC) && ...
                           length(tem_info.u_DC) == n_dof);
                    u(:, tt) = u(:, tt) + M*tem_info.u_DC;
                end
                switch tem_info.src_info.type
                    % -> Note that we use the same function handle b_fun to
                    %    represent the amplitude of either j or dj/dt!
                    % -> Note that the RHS vector b is associated with the
                    %    spatial distribution / geometric representation
                    case 'shut-off'
                    % Get initial value, i.e. current density j.
                    % -> M*u_0 = j(t=0) + ...
                    u(:, tt) = u(:, tt) + b*b_fun(t(tt));

                    case {'shut-off-gauss', 'ramp'}
                    % Incorporate rhs function, i.e. time derivative if
                    % current density dj/dt.
                    % -> (M + dt*K)*u_1 = ... + dt*dj/dt
                    u(:, tt) = u(:, tt) + dt(tt)*b*b_fun(t(tt));

                end
                % -> u_1 = (M + dt*K)^-1 * [M*u0 + dt*dj/dt]
                u(:, tt) = solver_MdtK{solver_idx_in_t(tt)}.solve(u(:, tt));
            else
                % Forward propagation in time:
                % -> u_i = (M + dt*K)^-1 * [M*u_i-1 + dt*dj/dt]
                u(:, tt) = solver_MdtK{solver_idx_in_t(tt)}.solve(M*u(:, tt-1) + ...
                                                            dt(tt)*b*b_fun(t(tt)));
            end
            fprintf(repmat('\b', 1, str_verbose));

            % Transform u to observation quantity (only keep t_obs values).
            if t_obs_idx_in_t(tt) ~= 0
                % -> d_i^obs = O * u_i
                d_(:, t_obs_idx_in_t(tt)) = (O.' * u(:, tt)).';
            end
        end
        fprintf('%.2ds \n', toc(run_time));

        % Reshape into vertical concatenated expression, i.e. the data vector.
        d = d_(:);

        %%% Calculate sensitivity via backard propagation %%%

        if nargout > 1
            % Get sensitivity.
            if n_d >= n_dof/2 % arbitrary treshold!
                % Skip calculation of S would be too large.
                warning(['It seems that the sensitivity should be calculated ', ...
                         'for the FE coefficients which is extremly inefficient. ', ...
                         'Skip calculation.']);
                S = [];

            else
                % Prepare calculation.
                S = zeros(n_d, n_t_obs, n_param);

                % Calculate sensitivity.
                if tem_info.solver_parallel > 1
                    S = get_S_parallel(S);
                else
                    S = get_S_serial(S);
                end

                % Reshape into vertical concatenated expression (see d).
                S = reshape(S, [], n_param);
            end
        end

        %%% Clean up %%%

        % Remove all solver objects.
        for dd = 1:n_dt
            solver_MdtK{dd}.delete();
        end

        %% Helper.

        function b = get_djdt_fun(src_info)
            % Define function handle for time-dependet rhs (dj(t)/dt).

            % Fetch info
            switch src_info.type
                case {'ramp', 'shut-off-gauss'}
                    assert(isfield(src_info, 'param'));
            end

            % Define source function handle w.r.t. time.
            switch src_info.type
                case 'ramp'
                % Define square-wave impulse for 0 >= t <= ramp end ( = d[linear j]/dt).
                b = @(t) wrap_ramp(t, src_info.param);

                case 'shut-off'
                % Define dirac-shut-off at t=0 ( = d[Heaviside j]/dt).
                b = @(t) wrap_shut_off(t);

                case 'shut-off-gauss'
                % Define Gauss-impulse centered at t=0 (= finite length approx. to d[Heaviside j]/dt).
                b = @(t) 1/(src_info.param * sqrt(2*pi)) .* ...
                            exp(-1/2*((t)/src_info.param).^2);

                otherwise
                    error('Unknown src time derivative function');
            end

            % Helper.
            function b = wrap_ramp(t, src_param)
                if t < 0 || t > src_param
                    b = 0;
                else
                    b = 1/src_param; % assume unit src current
                end
            end

            function b = wrap_shut_off(t)
                if t == 0
                    b = 1;
                else
                    b = 0;
                end
            end
        end

        function [t, dt] = adjust_t_dt(tem_info)
            % Adjust t to numerically handle the source turn-off-function.
% FIXME: may better split up times t=[0, t0] (or [0, t_obs(1)]) into a
%        fixed number of steps

            % Fetch info.
            t0 = tem_info.t(1);
            t_ = tem_info.t;
            dt_ = tem_info.dt;
            dt_ref = tem_info.src_info.dt_ref;

            % Get number of points to add for dt_ref such that b_fun is approx. 0.
            n_t_add_ref = add_time_steps(dt_ref, b_fun);

            % Check that b_fun = 0 (including thresold steps) is reached before t0.
            % -> Ensure to have add_ref_tresh tiny time steps after initial src signal
            add_ref_tresh = 1e1;
            if (t0 - (dt_ref*(n_t_add_ref))) < 0
                error('dt_ref is too large, no time steps can be added at all.');
            elseif (t0 - (dt_ref*(n_t_add_ref + add_ref_tresh))) < 0
                warning(['Could not add %d additional steps between t0 and t=0 for given dt_ref. ', ...
                         'Proceed with 0 additional steps'], add_ref_tresh);
                add_ref_tresh = 0;
            else
                t_tmp = ceil((t0 - (dt_ref*n_t_add_ref)) / dt_ref);
                if add_ref_tresh > t_tmp
                    warning('Need to reduce to %d additional steps between t0 and t=0.', t_tmp);
                    add_ref_tresh = t_tmp;
                end
            end
            t_tmp = (n_t_add_ref + add_ref_tresh) * dt_ref;
            assert((t0-t_tmp) > 0);

            % Fill up time steps with original dt_ to reduce number of tiny time
            % steps to be added.
            n_t_add_dt = floor((t0-t_tmp) / dt_(1));
            t_ = [fliplr(cumsum([t0, zeros(1, n_t_add_dt) - dt_(1)])), t_(2:end)];
            dt_ = [zeros(1, n_t_add_dt) + dt_(1), dt_];

            % Slightly adjust dt_ref such that final time series step lengths only
            % consists of dt_ref and elements of dt_ (i.e. avoid to introduce new
            % dt at transition between dt_ref and dt_(1)).
            t_ref_steps = ceil(t_(1)/min(dt_ref, dt_(1))); % ensure to refer to smallest dt
            assert(t_ref_steps >= (n_t_add_ref + add_ref_tresh));
            % Add time steps in + dir.
            t_ref_pos = linspace(0, t_(1), t_ref_steps + 1);
            if n_t_add_dt == 0
                assert(t_ref_pos(end) == t0);
                t_ref_pos(end) = []; % make sure not to add duplicates
            end
            dt_ref = t_ref_pos(2) - t_ref_pos(1);

            % For perfect shut-off, only one point will be added and no neg. time
            % axis is required.
            if n_t_add_ref > 1
                % Get number of time steps added in -t dir. s.t. abs(db) < t_tol.
                n_t_add_ref = add_time_steps(dt_ref, b_fun);
                t_ref_neg = fliplr(cumsum([0, zeros(1, n_t_add_ref) - dt_ref]));
            else
                t_ref_neg = 0;
            end

            % Redefine t_ an dt_.
            t = [t_ref_neg, t_ref_pos(2:end), t_];
            dt = [dt_ref + zeros(1, length(t_ref_neg) + length(t_ref_pos)-1), dt_];
            assert(length(t) == length(dt));

            % Helper.
            function n_t_add = add_time_steps(dt, b)
                % Get number of time steps to be added for reduction of b(t) < tol.
                %
                % Assumptions:
                %   (i)  b(t) -> 0 for t -> +- inf
                %   (ii) signal starts at t = 0

                n_t_add = 1;
                tol = 3e-2;     % limit up to which b is decreased
                treshold = 250; % max. number of additional points
                while abs(b(n_t_add*dt)) > tol
                    n_t_add = n_t_add+ 1;
                    if n_t_add > treshold
                        warning('Too much additional steps, inrease tol.');
                        break;
                    end
                end
            end
        end

        function S = get_S_serial(S)

            % Loop over 'solver chunks', i.e. unique solver objects.
            fprintf('Backpropagate sensitivity for %d times in serial w.r.t. %d RHS: ', n_t_obs, n_d);
            run_time = tic;
            n_FBS = 0;
            for cc = 1:n_dt

                % Assign observation times t_obs belonging to current solver.
                t_obs_cur_idx_in_solver_and_t = solver_idx_in_t == cc & t_obs_idx_in_t;

                % Assign index of time steps t belonging to current t_obs.
                idx_t_obs_cur_in_t = find(t_obs_cur_idx_in_solver_and_t);
                n_t_obs_cur = length(idx_t_obs_cur_in_t);
                if n_t_obs_cur == 0
                    % Skip current chunk as no t_obs belonging.
                    continue;
                end

                % Loop over t_obs that share the same initial solver,
                % starting from latest time.
                if n_t_obs_cur >= 1
                    % Only if current solver chunk belongs to initial solver of
                    % more than one t_obs, keep adjoint solutions in memory.
                    idx_store_max = idx_t_obs_cur_in_t(end)-idx_t_obs_cur_in_t(1)+1;
                    U_adjoint_ = zeros(n_d, n_dof, idx_store_max);
                end
                for ii = n_t_obs_cur:-1:1

                    % Get adjoint solution via back propagation.
                    U_adjoint_rhs = O; % initial step refers to observation operator
                    S0 = zeros(n_d, n_param);
                    str_verbose = fprintf('%d', t_obs_idx_in_t(idx_t_obs_cur_in_t(ii)));
                    for tt_ = idx_t_obs_cur_in_t(ii):-1:1

                        % For the n_t_obs-1 t_obs in current chunk:
                        idx_U_adjoint = idx_store_max-(idx_t_obs_cur_in_t(ii) - tt_);
                        if (ii ~= n_t_obs_cur) && (tt_ >= idx_t_obs_cur_in_t(1))
                            % Load adjoint solutions.
                            U_adjoint = U_adjoint_(:, :, idx_U_adjoint);
                        else
                            % Calculate adjoint solution.
                            % -> calculation is also required for all times
                            %    that belong to a different solver chunk
                            %    than the current chunk in the outermost iteration
                            U_adjoint = solver_MdtK{solver_idx_in_t(tt_)}.solve(U_adjoint_rhs).';
                            n_FBS = n_FBS+1;

                           % For latest t_obs in current chunk:
                            if (ii == n_t_obs_cur) && (tt_ >= idx_t_obs_cur_in_t(1)) && ...
                               (n_t_obs_cur > 1)
                                % Store all adjoint solutions up to the time step,
                                % where dt (i.e. the solver) changes.
                                U_adjoint_(:, :, idx_U_adjoint) = U_adjoint;
                            end
                        end
                        if (tt_ <= idx_t_obs_cur_in_t(1)) || ...
                           (ii == n_t_obs_cur) && (tt_ >= idx_t_obs_cur_in_t(1))
                            % Update rhs.
                            % -> symmetric M does not change with time
                            U_adjoint_rhs = (U_adjoint*M).';
                        end

                        % Evaluate sensitivity.
                        if tt_ == 1
                            % -> u: (u_DC) - u_i
                            if isfield(tem_info, 'u_DC')
                                u_ = tem_info.u_DC - u(:, tt_);
                            else
                                u_ = -u(:, tt_);
                            end
                            S0 = (U_adjoint * ttv(TM, u_, 2)) + S0;
                        else
                            % -> u: u_i-1 - u_i
                            S0 = (U_adjoint * ttv(TM, u(:, tt_-1) - u(:, tt_), 2)) + S0;
                        end
                    end
                    S(:, t_obs_idx_in_t(idx_t_obs_cur_in_t(ii)), :) = S0;
                    fprintf(repmat('\b', 1, str_verbose));
                end
            end
            fprintf('%.2ds (%i FBS)\n', toc(run_time), n_FBS);
        end

        function S = get_S_parallel(S)

            % Destroy (local) mumps instances.
            for i = 1:n_dt
                switch solver_type
                    case 'mumps_ooc'
                        % Keep meta info but destroy Mumps object
                        solver_MdtK{i}.unload();
                    case 'chol'
                        % Ensure that no handle will be serialized.
                        solver_MdtK{i} = rmfield(solver_MdtK{i}, 'solve');
                    otherwise
                        solver_MdtK{i}.delete();
                end
            end

            % Define pool constants for all large, read-only variables too,
            % such that they are copied only once.
            worker_solver_idx_in_t = parallel.pool.Constant(solver_idx_in_t);
            worker_u = parallel.pool.Constant(u);
            worker_TM = parallel.pool.Constant(sol.TM);
            if isfield(tem_info, 'u_DC')
                worker_u_DC = parallel.pool.Constant(tem_info.u_DC);
            else
                worker_u_DC = [];
            end
            worker_M = parallel.pool.Constant(M);
            worker_K = parallel.pool.Constant(K);
            worker_O = parallel.pool.Constant(O);
            % Var1: (get_local_solver needs to be adjusted, respectively!)
            % - seems to be little bit slower than Var2
            % worker_solver_MdtK = parallel.pool.Constant(@() get_local_solver(solver_MdtK, K, M));
            % Var2: (get_local_solver is correctly set using .Value)
            worker_solver_MdtK = parallel.pool.Constant(@() get_local_solver(solver_MdtK));

            % Set current back-propagation time.
            t_obs_idx = find(t_obs_idx_in_t);

            % Loop over t_obs in parallel.
            fprintf('Backpropagate sensitivity for %d times in parallel w.r.t. %d RHS: ', n_t_obs, n_d);
            run_time = tic;
            parfor (ii = 1:n_t_obs, parpool_opts)
            % parfor ii = 1:n_t_obs
                % Set initial step which refers to observation operator.
                all_local_solver = worker_solver_MdtK.Value;
                solver2t_map = worker_solver_idx_in_t.Value;
                local_TM = worker_TM.Value;
                local_M = worker_M.Value;
                local_u = worker_u.Value;
                local_O = worker_O.Value;
                if isempty(worker_u_DC)
                    local_u_DC = 0;
                else
                    local_u_DC = worker_u_DC.Value;
                end
                S0 = zeros(n_d, n_param);
                % Get adjoint solution via back propagation.
                for tt_ = t_obs_idx(ii):-1:1
                    % Get adjoint solution.
                    % -> Pick the correct solver for this time step
                    local_solver = all_local_solver{solver2t_map(tt_)};
                    U_adjoint = local_solver.solve(local_O).';

                    % -> symmetric M does not change with time
                    local_O = (U_adjoint * local_M).';

                    % Evaluate sensitivity.
                    if tt_ == 1
                        % -> u: (u_DC)-u_i
                        u_ = local_u_DC - local_u(:, tt_);
                        S0 = (U_adjoint * ttv(local_TM, u_, 2)) + S0;
                    else
                        % -> u: u_i-1 - u_i
                        S0 = (U_adjoint * ttv(local_TM, local_u(:, tt_-1) - ...
                                              local_u(:, tt_), 2)) + S0;
                    end
                end
                S(:, ii, :) = S0;
            end
            fprintf('%.2ds (%i FBS)\n', toc(run_time), sum(t_obs_idx));

            % Helper
            function local_solver_MdtK = get_local_solver(local_solver_MdtK)
                % Load / create local factorizations (only once at startup)

                for jj = 1:n_dt
                    switch solver_type
                        case 'mumps_ooc'
                        local_solver_MdtK{jj}.load();
                        otherwise
                        local_solver_MdtK{jj} = get_solver_fun(dt_unique(jj)*worker_K.Value + worker_M.Value, solver_type);
                    end
                end
            end
        end
    end
end
