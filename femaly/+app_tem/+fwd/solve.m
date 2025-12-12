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
    %
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
                solver.delete = @() clear('solver_'); %required?

            case  {'mumps', 'mumps_ooc'}
                id_.SYM = 0;
                if issymmetric(A)
% FIXME: Why symmetry argument leads to seg.fault error?
                    % id_.SYM = 1;
                end
                if strcmp(type, 'mumps')
                    solver = solving.MUMPS(A, id_);
                else
                    % Choose solver variant.
                    solver = struct();
                    if isreal(A)
                        solver.mumps = @dmumps;
                    else
                        solver.mumps = @zmumps;
                    end

                    % Set file name var env.
                    assert(isscalar(num));
                    solver.num = num2str(num);
                    setenv('MUMPS_SAVE_PREFIX', ['deco_', solver.num]);

                    % Initialize mumps object.
                    solver.id = dmumps(initmumps());
                    solver.id.SYM = id_.SYM;
                    
                    % Set parameter to enable OOC.
                    solver.id.ICNTL(22) = 1;
                    solver.id.ICNTL(35) = 3;
                    
                    % Factorize OOC.
                    solver.id.JOB = 4;
                    solver.id = dmumps(solver.id, A);
                    
                    % Save the instance to disk.
                    solver.id.JOB = 7;
                    solver.id = dmumps(solver.id, A);
                    n_sys = size(A, 1);
                    solver.A = sparse([], [], [], n_sys, n_sys); % dummy

                    % Free internal memory (saved files remain).
                    solver.id.JOB = -2;
                    [~] = dmumps(solver.id, solver.A);

                    % Set wrapper.
                    solver.solve = @(b) solve_ooc(solver, b);
                    solver.delete = @() solve_delete(solver); 
                end

            otherwise
                error('Unsupported solver type.');
        end

        function x = solve_ooc(solver, b)
            % Initialize.
            id = dmumps(initmumps());
            setenv('MUMPS_SAVE_PREFIX', ['deco_', solver.num]);
            
            % % Reset solver parameter.
            id.SYM = solver.id.SYM;
            
            % Restore saved instance.
            id.JOB = 8;
            id = solver.mumps(id, solver.A);

            % Set RHS and solve.
            id.RHS = b;
            id.JOB = 3;
            id = solver.mumps(id, solver.A);
            x = id.SOL;

            % Cleanup object.
            id.JOB = -2;
            [~] = solver.mumps(id, solver.A);
        end

        function solve_delete(solver)
            % Initialize.
            id = solver.mumps(initmumps());
            setenv('MUMPS_SAVE_PREFIX', ['deco_', solver.num]);
            
            % Restore saved instance.
            id.JOB = 8;
            id = solver.mumps(id, solver.A);
            
            % Remove saved files.
            id.JOB = -3;
            id = solver.mumps(id, solver.A);
            
            % Cleanup object.
            id.JOB = -2;
            [~] = solver.mumps(id, solver.A);
        end
    end

    function [d, S] = solve_rba(b, K, M, O, solver_type, tem_info)
        % Create forward solution and sensitivity handle.

        % Get handle to forward solver.
        if tem_info.solver_parallel > 1
            warning('Not implemented yet.');
        end
        solver_fun_ = @(A) get_solver_fun(A, solver_type);

        % Initialize solver for mass matrix.
        solver_fun = solver_fun_(M);

        % Get initial value.
        x0 = solver_fun.solve(b);
        solver_fun.delete;

        % Solve.
        t = tem_info.t_obs;
        d = zeros(size(O, 2), length(t));
        u = zeros(size(M, 1), length(t));
        r = rba.create_rba_matrix(tem_info.n_poles, K, M, x0, @(A, b) wrap_solver(A, solver_fun_, b));
        fprintf('Calculate RBA for %d time steps and %d DOF: ', length(t), size(M, 1));
        run_time = tic;
        for tt = 1:length(t)
            % Calculate time solution.
            str_verbose = fprintf('%d', tt);
            u(:, tt) = r(t(tt));
            d(:, tt) = O.' * u(:, tt);
            fprintf(repmat('\b', 1, str_verbose));
        end
        fprintf('%.2ds \n', toc(run_time));

        % Clean up.
        solver_fun.delete;

        % Reshape into vertical concatenated expression.
        d = d(:);

% TODO: implement from Marios Habil.
        if nargout > 1
            S = [];
            warning('Not implemented yet.');
        end

        %%% Helper %%%

        function u = wrap_solver(A, solver_fun_, b)
            % Solve fwp and crude overwrite output.

            assert(isvector(b) && size(b, 1) == size(A, 2));
            solver = solver_fun_(A);
            u = solver.solve(b);
            solver.delete;
        end
    end

    function [d, S] = solve_impl(b, K, M, O, TM, solver_type, tem_info)
        % Create forward solution and sensitivity matrix.
        %
        % t = 0 is assumed as initial time
        
        % Get forward problem solver.
        if tem_info.solver_parallel > 1
            if ~contains(solver_type, 'mumps')
                warning(['Parallel computing requires mumps solver. ', ...
                         'Proceed with MUMPS_ooc.'])
            end
% HACK: handle cases
solver_type = 'mumps_ooc';

            % Set directories & prefixes.
            save_dir = fullfile(pwd, 'tmp_mumps');  % storage location
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);    
            end
            
            % Set required env vars.
            if isempty(getenv('MUMPS_SAVE_DIR'))
                setenv('MUMPS_SAVE_DIR', save_dir);
            end
            solver_fun = @(A, num) get_solver_fun(A, solver_type, num);
        else
% HACK: handle cases
if contains(solver_type, 'mumps')
solver_type = 'mumps';
end
            solver_fun = @(A, num) get_solver_fun(A, solver_type);
        end

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

        % Get all unique solver objects.
        solver_MdtK = cell(n_dt, 1);
        fprintf('Calculate %d system decompositions for %d DOF: ', n_dt, size(M, 1));
        run_time = tic;
        for dd = 1:n_dt
            str_verbose = fprintf('%d', dd);
            solver_MdtK{dd} = solver_fun(dt_unique(dd)*K + M, dd);
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
        % -> only for sensitivity calculation parallelism can be exploited
            if ~isfield(tem_info, 'solver_parallel')
                tem_info.solver_parallel = 0;
            elseif tem_info.solver_parallel == 1
                warning(['Only 1 worker should be used. ', ...
                        'Proceed with serial calculation.']);
                tem_info.solver_parallel = 0;
            elseif tem_info.solver_parallel > 1
                % Check if parallel toolbox license is available.
                if license('test', 'distrib_computing_toolbox')
                    % Check for existing pool of parallel worker.
                    par_pool = gcp('nocreate');
                    if isempty(par_pool)
                        % Create pool.
                        parpool('local', tem_info.solver_parallel);
                    else
                        fprintf(['Using active pool with %d worker. ', ...
                                'Use delete(gcp(''nocreate'')); to close current pool.\n'], par_pool.NumWorkers);
                    end
                else
                    warning(['No Distributed Computing Toolbox available. ', ...
                            'Proceed with serial calculation.']);
                    tem_info.solver_parallel = 0;
                end
            end
            
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
                    S = get_S_parallel_ooc(S);
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
            if isa(solver_MdtK, 'parallel.pool.Constant')
                solver_MdtK.Value{dd}.delete;
            elseif isa(solver_MdtK, 'cell')
                solver_MdtK{dd}.delete();
            else
                error('Unknown type.');
            end
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

            % Fetch info.
            t0 = tem_info.t(1);
            t_ = tem_info.t;
            dt_ = tem_info.dt;
            dt_ref = tem_info.src_info.dt_ref;

            % Get number of points to add for dt_ref such that db is approx. 0.
            n_t_add_ref = add_time_steps(dt_ref, b_fun);

            % Check that db = 0 (including thresold steps) is reached before t0.
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

        function S = get_S_parallel_ooc(S)
                
            % Define pool constants for all large, read-only variables too, 
            % such that they are copied only once.
            worker_u = parallel.pool.Constant(u);
            worker_solver_idx_in_t = parallel.pool.Constant(solver_idx_in_t);
            worker_TM = parallel.pool.Constant(sol.TM);
            if isfield(tem_info, 'u_DC')
                worker_u_DC = parallel.pool.Constant(tem_info.u_DC);
            else
                worker_u_DC = [];
            end
            worker_M = parallel.pool.Constant(M);
            worker_O = parallel.pool.Constant(O);
            worker_solver_MdtK = parallel.pool.Constant(solver_MdtK);

            % Set current back-propagation time.
            t_obs_idx = find(t_obs_idx_in_t);

            % Loop over t_obs in parallel.
            fprintf('Backpropagate sensitivity for %d times in parallel w.r.t. %d RHS: ', n_t_obs, n_d);
            run_time = tic;
            parfor ii = 1:n_t_obs
                % Set initial step which refers to observation operator.
                U_adjoint_rhs = worker_O.Value; 
                S0 = zeros(n_d, n_param);
                % Get adjoint solution via back propagation.
                for tt_ = t_obs_idx(ii):-1:1
                    % Get adjoint solution.
                    % -> Pick the correct solver for this time step
                    U_adjoint = worker_solver_MdtK.Value{worker_solver_idx_in_t.Value(tt_)}.solve(U_adjoint_rhs).';

                    % -> symmetric M does not change with time
                    U_adjoint_rhs = (U_adjoint * worker_M.Value).';

                    % Evaluate sensitivity.
                    if tt_ == 1
                        % -> u: (u_DC)-u_i
                        if ~isempty(worker_u_DC)
                            u_ = worker_u_DC.Value - worker_u.Value(:, tt_);
                        else
                            u_ = -worker_u.Value(:, tt_);
                        end
                        S0 = (U_adjoint * ttv(worker_TM.Value, u_, 2)) + S0;
                    else
                        % -> u: u_i-1 - u_i
                        S0 = (U_adjoint * ttv(worker_TM.Value, worker_u.Value(:, tt_-1) - ...
                                              worker_u.Value(:, tt_), 2)) + S0;
                    end
                end
                S(:, ii, :) = S0;
            end
            fprintf('%.2ds (%i FBS)\n', toc(run_time), sum(t_obs_idx));
        end
    end
end
