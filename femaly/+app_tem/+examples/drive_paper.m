% TD-EM problem for a rectangular loop at top of hom. halfspace
%
% PDE:
%   curl (1/mu curl e) + sigma d_t(e) = d_j/d_t (e = e(x, t))
% where:
%   d_j/d_t ~= 0 only if a source ramp time is considered
%   IC: sigma e(t=0)  = j
% where:
%             j      ~= 0   if   d_j/d_t = 0
%   BC:       n x e   = 0       (n = surface normal)
%
% Weak form (involving Stokes)
%   int_V 1/mu curl(Phi_i) * curl(Phi_j) dx
%     + d_t int_V sigma  Phi_i * Phi_j dx = int_V d_j/d_t * Phi_j dx
%
% IC:
%   int_V sigma  Phi_i * Phi_j dx = int_V j * Phi_j dx
%
% lin. system
%   Ku + Mu = b * f(t)
% IC:
%   Mu = b
% with
%   b = sum_e p_e * int_L Phi_i * tau ds * f(t)
% and
%   p_e = I dl dipole moment for current I on length dl
%   tau        edge tangent
%   f(t) =     a time-depending amplitude modulation
%
% Implicit Euler:
%   -M*u(ti) + (M + dt*K)u(ti+1) = b
%                        u(ti+1) = (M + dt*K)^-1 * (M*u(ti) + b * f(t))
%                        u(t1)   = (M + dt*K)^-1 * (b       + b * f(t1))

% Jacobian (Frechet derivatives) for implicit Euler:
%   d_u(t_i+1)/d_m = (M + dt*K)^-1 [d_M/d_m*(u_i - u_i+1) + M*d_u_i/d_t]
%   d_u(t_1)/d_m   = (M + dt*K)^-1 [d_M/d_m*(-u_1)]               as   d_b/d_m = 0

%% Initialize.

% Start script from within femaly main folder!

export_results = util.pick(1, false, true);

% Define equation solver.
% Note: only for implicit Euler, a time-depending source function can be
%       defined.
integration_type = util.pick(1, 'impl');
integration_times = util.pick(2, 'equal', 'split'); % only acts on implicit euler
solver_type = util.pick(1, 'mumps_ooc', 'mumps', 'chol');
solver_parallel = 0;                                % number of matlab worker (>1 for parallel computing)
src_type = util.pick(1, 'shut-off', 'ramp');
data_type = util.pick(3, {'dBdt'}, {'E'}, {'E', 'dBdt'});

% Define FE parameter.
fe_order = 2;

% Verbosity parameter.
file_path = ['results', filesep, 'tem_paper', filesep];
% diary([file_path, 'data.txt']);

% Domain parameter
sig_air = 1e-9;
sig_earth = 1/10;
%
str_param = sprintf(' sig_air: %.2d\n sig_earth: %.2d\n mue_0: %d\n', ...
                      sig_air, sig_earth, util.EMConstants.mu_0);

% Meshing parameter.
shape = util.pick(1, 'sphere', 'cuboid');
domain_r = 1e5;
src_length = util.pick(1, 5, 1);
size_at_wr = src_length/10;
size_at_pt = 5e-1;
size_at_box = util.pick(2, 1e0, domain_r);
ref_global = 0;
line = [-src_length/2, -src_length/2, 0, 1, 1;    % square, centered at [0,0,0]
         src_length/2, -src_length/2, 0, 1, 1;
         src_length/2,  src_length/2, 0, 1, 1;
        -src_length/2,  src_length/2, 0, 1, 1;
        -src_length/2, -src_length/2, 0, 1, 1, ...
        ];
src_A = 1/2 * sum(line(1:end-1, 1).*line(2:end, 2) - line(2:end, 1).*line(1:end-1, 2));

point = [20, 0, 0, 1, 1];       % setup Strokkur

% Time / solver paremeter.
t_exp_min = log10(1e-7);
t_exp_max = log10(1e-3);
t0 = 10^t_exp_min;

% Define gate times t_obs_.
%{
t_sys_ = [1.81e-6, 3.69e-6, 6.1e-6, 9.2e-6, 1.31e-5, 1.81e-5, 2.44e-5, ...
         3.26e-5, 4.29e-5, 5.61e-5, 7.29e-5, 9.44e-5, 1.217e-4, 1.614e-4, ...
         2.01e-4, 2.578e-4, 3.302e-4, 4.227e-4, 5.406e-4, 6.909e-4]; % time after t_0
t_0_ = 0.9e-6;                                                       % time where TX ramp is assumed to be 0
t_obs_ = t_sys_ + t_0_;                                               % approximately values from measurement
%}
t_obs_ = logspace(t_exp_min, t_exp_max, 20);                          % approximately logarthmically equidistant values from RBA Setting

% Set time integration parameter
switch integration_type
    case 'rba'

    % Set up observations times.
    n_poles = 12;
    %
    str_int = sprintf(' t_min: exp(%.1f), t_max: exp(%.1f)\n number of time steps: %d\n poles: %d\n', ...
                      t_exp_min, t_exp_max, t_steps, n_poles);

    t = t_obs; %Dummy - here calculated times are exactly the observed ones

    case 'impl'
    % Set general time vector.
    switch integration_times
        case 'equal'
            t_steps = 175;          % per decade
            [t, dt] = equal_spaced_decades(t_exp_min, t_exp_max, t_steps);
        case 'split'
            %{
            t_steps = 7;            % per splitting intervall
            t_num_split = 17;
            [t, dt] = splitted_intervalls(logspace(t_exp_min, t_exp_max, t_num_split), t_steps);
            %}
            t_steps = 27;
            [t, dt] = splitted_intervalls(t_obs_, t_steps);
        otherwise
            error('Unknown time axis definition.');
    end

    % Define time-derivative of source function, i.e. dj/dt.
    switch src_type
        case 'shut-off'
        % Dirac         -> derivative Heavyside function
        db_param = [];
        dt_ref = t0/2e1;
        %
        str_src = sprintf(' src_type: %s\n dt_ref at 0: %.2d\n', src_type, dt_ref);

        case 'ramp'
        % Square-wave   -> derivative of perfect linear ramp
        ref_portion = 10;
        db_param = 5e-8;                    % ramp length (Geonics, 10x10 loop)
        dt_ref = db_param/ref_portion;
        %
        str_src = sprintf(' src_type: %s\n dt_ref for ramp: %.2d\n ramp length: %.2d\n', ...
                           src_type, dt_ref, db_param);
    end
    %
    str_int = [sprintf(' t0: exp(%.1f), t_max: exp(%.1f)\n dt(1): %.2d  -> steps/decade: %d\n', ...
                       t_exp_min, t_exp_max, dt(1), t_steps), str_src];

    % Find values in t closest to t_obs.
    % -> if t_split ~= t_obs_ values won't be simulated exactly
    [~, idx] = arrayfun(@(x) min(abs(t - x)), t_obs_);

    % Select unique (in case some values double) the corresponding times.
    t_obs = unique(t(idx));
    t = t(1:idx(end));                      % throw away times > max(t_obs)
    dt = dt(1:idx(end));
end
%
str_tmp = arrayfun(@(i) {[sprintf('%d ', point(i, 1:3)), ' | ']}, 1:size(point, 1));
str_geo = sprintf(' shape: %s\n domain_r: %d\n src_length (center: [0,0,0]): %d\n src_A: %d\n point: %s\n', ...
                   shape, domain_r, src_length, src_A, [str_tmp{:}]);
str_meshing = sprintf(' size_at_wr: %.2d\n size_at_pt: %.2d\n', ...
                      size_at_wr, size_at_pt);

%% Get mesh and parameter vector.

% Run meshing.
keep_files = true;
keep_air = true;
[mesh, cm, fm, lm, pm, phy] = meshing.generate_mesh3D(...
                        'domain_r', domain_r, ...
                        'line', line, ...
                        'point', point, ...
                        'shape', shape, ...
                        'keep_air', keep_air, ...
                        'keep_files', keep_files, ...
                        'size_at_wr', size_at_wr, ...
                        'size_at_pt', size_at_pt, ...
                        'ref', ref_global, ...
                        'geo_code', meshing_code([0, 20], [-5, 5], [-5, 0], 10, size_at_box), ...
                        'marker', [3, 2, 1, 0, -1]);
mesh.init_geometric_queries;
mesh.compute_boundary_facets;
c_mp = mesh.get_cell_centroids;
str_mesh = sprintf(' n_cell: %d\n', size(mesh.cells, 2));

% Identify boundary and line.
get_id = @(ent, name) phy{ent}{1}(ismember(phy{ent}{2}, name)).';
line_m = get_id(3, 'line_1'); % or use id defined above for identification
air_m = get_id(1, 'air');
earth_m = get_id(1, 'halfspace');

% Set parameter vector.
sigma = zeros(size(mesh.cells, 2), 1);
sigma(cm == air_m) = sig_air;
sigma(cm == earth_m) = sig_earth;
assert(~(any(sigma == 0)));

%% Prepare system assembly.

% Prepare BC info (hom. Dirichlet everywhere).
BC_fun_ = {{@(x) true}, {@(x) zeros(1, mesh.dim)}};
bnd_sum = {'dirichlet', BC_fun_};
%
tem_info = struct();
if strcmp(integration_type, 'impl')
    src_info = struct();
    src_info.type = src_type;
    src_info.dt_ref = dt_ref;
    src_info.param = db_param;
    tem_info.src_info = src_info;
    tem_info.dt = dt;
else
    tem_info.n_poles = n_poles;
end
point_obs = point(1, :);                          % just consider first of all points
tem_info.pt_RX = point_obs;
tem_info.pt_TX = line;
tem_info.t_obs = t_obs;
tem_info.t = t;
tem_info.lm = lm;
tem_info.srcm = line_m;
tem_info.integration_type = integration_type;
tem_info.integration_times = integration_times;
tem_info.solver_parallel = solver_parallel;
tem_info.data_type = data_type;

% Assemble ans solve.
sol = app_tem.fwd.assemble('3D', mesh, fe_order, bnd_sum, tem_info);
str_fe = sprintf(' n_dof: %d\n fe_order: %d\n', sol.dofmap.dim, fe_order);
%
wrap_solver = @(sigma) app_tem.fwd.solve(sol, sigma, solver_type, tem_info);
%
timing = tic;
[dB, S] = wrap_solver(sigma);
tmp = toc(timing);
str_timing = sprintf(' runtime: %.2f\n', tmp);
[~, line_idx] = find(t_obs(:) == t(:).');
n_d = size(sol.O, 2);

% Check S.
if sol.dofmap.dim < 1e4 && ~isempty(S)
    figure(1);
    clf;
    util.run_taylor_test(wrap_solver, sigma, S);
end

fprintf([sprintf('Integration (%s):\n', integration_type), str_int, '\n']);
fprintf(['Geometry: \n', str_geo, '\n']);
fprintf(['Parameter:\n', str_param, '\n']);
fprintf(['Mesh:\n', str_meshing, str_mesh, '\n']);
fprintf(['FE: \n', str_fe, str_timing, '\n'])
% diary off

%% Compare and visualize.

% Loop over observations.
fig_num = 2;
dB_fe = reshape(dB, mesh.dim, size(tem_info.pt_RX, 1), length(tem_info.data_type), length(tem_info.t_obs));
for dd = 1:length(tem_info.data_type)
    for pp = 1:size(tem_info.pt_RX, 1)
        assert(~isempty(point(pp, :)));

        % Calculate (semi)analytical reference solutions at !original
        % t_obs!.
        %
        % Get reference solution from empymod (finite-length wire TX).
        scale = util.EMConstants.mu_0;
        mrec = 'dH';
        if ~isfield(tem_info, 'data_type') || strcmp(tem_info.data_type{dd}, 'dBdt')
            % Nothing to do
            vmd_type = {'dHr', 'dHz'};
        elseif strcmp(tem_info.data_type{dd}, 'dHdt')
            scale = 1;
            vmd_type = {'Hr', 'Hz'};
        elseif strcmp(tem_info.data_type{dd}, 'E')
            vmd_type = {'Ef'};
            mrec = 'E';
            scale = 1;
        end

        % dH_wi = -reference.EH4wire(line(:, 1:3), [point(pp, 1:2), point(pp, 3)], [1/1, 1/100], [0, -0.6], ...
        dH_wi = -reference.EH4wire(line(:, 1:3), [point(pp, 1:2), point(pp, 3)], sig_earth, 0, ...
                                         'tx_type', 'E', 'rx_type', mrec, ...
                                         't_type', 'time', 't_list', t_obs_, ...
                                         'sigma_air', sig_air);
        dB_wi = dH_wi / src_A * scale;

        % Show observations - all scaled to unit strength/length.
        fig = figure(fig_num);
        clf;
        comp = {'x', 'y', 'z'};
        x_bnd = [min(t), max(t)];
        for cc = 1:3
            d_fe = (squeeze(dB_fe(cc, pp, dd, :))); % at t_obs closest to t_obs_
            d_emp = (squeeze(dB_wi(cc, :)));        % at t_obs_
            leg_str = {'fe-', 'fe+', 'empy-', 'empy+'};
            leg_use = [any(d_fe < 0), any(d_fe > 0), ...
                       any(d_emp < 0), any(d_emp > 0)];
            subplot(2, 3, cc)
                semilogx(t_obs(d_fe < 0), abs(d_fe(d_fe < 0)), 'r.', ...
                         t_obs(d_fe > 0), d_fe(d_fe > 0), 'b.', ...
                         t_obs_(d_emp < 0), abs(d_emp(d_emp < 0)), 'ro', ...
                         t_obs_(d_emp > 0), d_emp(d_emp > 0), 'bo');
                legend(leg_str{leg_use}, 'Location', 'SouthEast');
                if cc == 2
                    title({['type: ', integration_type, ...
                            sprintf('; RX: [%d, %d, %d]', point(pp, 1), point(pp, 2), point(pp, 3))], ...
                          [tem_info.data_type{dd}, '_', comp{cc}]});
                else
                    title([tem_info.data_type{dd}, '_', comp{cc}]);
                end
                xlabel('s');
                if cc == 1
                    ylabel('V/m²');
                end
                xlim(x_bnd);
                % ylim([1e-12, 1e-4]);
            subplot(2, 3, 3+cc)
                semilogx(t_obs(:), 100*(abs(d_fe(:))-abs(d_emp(:)))./abs(d_emp(:)), '.', ...
                         t_obs(:), t_obs(:)*0, '--');
                xlim(x_bnd);
                ylim([-20, 20]);
            if cc == 1
                ylabel('%');
            elseif cc == 2
                title('rel. err');
            end
        end
        if export_results
            savefig(fig, sprintf('plot_%i', fig_num));
        end
        fig_num = fig_num + 1;
    end
end

%% Export

if export_results

    % Calculate FE coefficients separately (for visualization of fields).
    n_vtx = size(sol.dofmap.mesh.vertex_coords, 2);
    n_t_obs = length(tem_info.t_obs);
    sol_ = sol;
    sol_.O = speye(sol.dofmap.dim); % HACK: calculate all dof for t_obs
    u = app_tem.fwd.solve(sol_, sigma, solver_type, tem_info);
    u = reshape(u, [], n_t_obs);

    % Export dB/dt for all mesh nodes.
    Q_all = assembling.assemble_point_sources(sol.dofmap, sol.dofmap.mesh.vertex_coords, 'curl');
    dB_dt = reshape(Q_all.' * u, 3, n_vtx, n_t_obs);
    file_str = 'dB_dt';
    for tt = 1:n_t_obs
        % Ensure leading 0 for sorting in dir(...)!
        io.writeVTU(sprintf([[file_path, file_str], '_%0', num2str(ceil(log10(max(1, abs(n_t_obs)+1)))), 'i.vtu'], tt), mesh, ...
                        'cell_marker', sigma, ...
                        'domain_marker', cm, ...
                        'point_marker', -dB_dt(:, :, tt));
    end
    % Update file list of current folder.
    tmp = dir(file_path);
    tmp = {tmp.name};
    io.writePVD([file_path, file_str], tmp(startsWith(tmp, file_str) & ~contains(tmp, {'.pvd', 'png', 'pdf', 'jpg', 'jpeg'})), t_obs);

    % Export E for all mesh nodes.
    Q_all = assembling.assemble_point_sources(sol.dofmap, sol.dofmap.mesh.vertex_coords);
    E = reshape(Q_all.' * u, 3, n_vtx, n_t_obs);
    file_str = 'E';
    for tt = 1:n_t_obs
        io.writeVTU(sprintf([[file_path, file_str], '_%0', num2str(ceil(log10(max(1, abs(n_t_obs)+1)))), 'i.vtu'], tt), mesh, ...
                        'cell_marker', sigma, ...
                        'domain_marker', cm, ...
                        'point_marker', -E(:, :, tt));
    end
    tmp = dir(file_path);
    tmp = {tmp.name};
    io.writePVD([file_path, file_str], tmp(startsWith(tmp, file_str) & ~contains(tmp, {'.pvd', 'png', 'pdf', 'jpg', 'jpeg'})), t_obs);

    % Export S.
    % -> weighted by cell volume to get physical correct representation
    n_c = length(sigma);
    n_dd = length(data_type);
    c_vol = mesh.get_cell_volumes;
    S_5D = reshape(S, mesh.dim, length(tem_info.data_type), size(tem_info.pt_RX, 1), length(tem_info.t_obs), n_c);
    % Export sensitivity for each vector component w.r.t. all observations.
    for dd = 1:n_dd
        for cc = 1:3
            file_str = sprintf('S_%s_%i', data_type{dd}, cc);
            for tt = 1:n_t_obs
                S_comp = squeeze(S_5D(cc, dd, :, tt, :)).' ./ c_vol.';
                io.writeVTU([[file_path, file_str], sprintf(['_%0', num2str(ceil(log10(max(1, abs(n_t_obs)+1)))), 'i.vtu'], tt)], mesh, ...
                                'cell_marker', S_comp, ...
                                'domain_marker', cm);
            end
            tmp = dir(file_path);
            tmp = {tmp.name};
            io.writePVD([file_path, file_str], tmp(startsWith(tmp, file_str) & ~contains(tmp, {'.pvd', 'png', 'pdf', 'jpg', 'jpeg'})), t_obs);
        end
    end

    % Export points.
    util.pt2txt(tem_info.pt_RX(:, 1:3), [file_path, 'pt_RX']);
    util.pt2txt(tem_info.pt_TX(:, 1:3), [file_path, 'pt_TX']);
end

%% Helper

function [t, dt] = equal_spaced_decades(pot_min, pot_max, sep)
    % Split up each decade of range 10^pot_min -  10^pot_max in sep pieces.

    % Sanity check.
    tmp = (floor(pot_min)==ceil(pot_min) && floor(pot_max)==ceil(pot_max));
    if ~tmp
        warning('Some boundary does not match 10^i, cutting off t at closest values.');
    end

    % Calculate times and time step sizes.
    pot = floor(pot_min):ceil(pot_max);
    [t, dt] = deal([]);
    for ii = 1:length(pot)-1
        dt_ = (10^pot(ii+1)-10^pot(ii)) / sep;
        dt = [dt, dt_]; %#ok<*AGROW>
        t = [t, 10^pot(ii):dt_:10^pot(ii+1)];
    end
    t = unique(t);

    % Bloat unique dt to vector of same size as t.
    dt = [reshape(repmat(dt, sep, 1), [], 1).', dt(end)];

    % Cut off arrays.
    if ~tmp
        tmp = t < 10^(pot_min) | t > 10^(pot_max);
        t(tmp) = [];
        dt(tmp) = [];
    end
end

function [t, dt] = splitted_intervalls(t_split, sep)
    % Define linear spacings with sep steps between elements of t_split.

    % Sanity check.
    assert(length(t_split) >= 2);
    assert(all(diff(t_split) > 0));
    assert(length(sep) == length(t_split)-1 || isscalar(sep));
    if isscalar(sep)
        sep = zeros(length(t_split) - 1, 1) + sep;
    end

    % Calculate times and time step sizes.
    [t, dt] = deal([]);
    for ii = 1:length(t_split)-1
        dt_ = (t_split(ii+1) - t_split(ii)) / sep(ii);
        dt = [dt, dt_]; %#ok<*AGROW>
        t = [t, t_split(ii):dt_:t_split(ii+1)];
    end
    t = unique(t);

    % Bloat unique dt to vector of same size as t.
    dt = [cell2mat(arrayfun(@(x, y) {zeros(x, 1) + y}, sep(:), dt(:))).', dt(end)];
end

function geo_code = meshing_code(xbo, ybo, zbo, tbo, size_at_box)
    % Define additional geo code for specific refinements.

    geo_code = [ ...
    'Field[6] = Box;', newline, ...
    ['Field[6].VIn = ', num2str(size_at_box), ';'], newline, ...
    'Field[6].VOut = domain_r/5;', newline, ...
    ['Field[6].XMin = ', num2str(xbo(1)), ';'], newline, ...
    ['Field[6].XMax = ', num2str(xbo(2)), ';'], newline, ...
    ['Field[6].YMin = ', num2str(ybo(1)), ';'], newline, ...
    ['Field[6].YMax = ', num2str(ybo(2)), ';'], newline, ...
    ['Field[6].ZMin = ', num2str(zbo(1)), ';'], newline, ...
    ['Field[6].ZMax = ', num2str(zbo(2)), ';'], newline, ...
    ['Field[6].Thickness = ', num2str(tbo),';'], newline, ...
    'Field[200] = Min;', newline, ...
    'Field[200].FieldsList = {10, 20, 6};', newline, ...
    'Background Field = {200};', newline, ...
    ];
end
