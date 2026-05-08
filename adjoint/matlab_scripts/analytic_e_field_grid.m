% Calculate the TEM electric field of a vertical magnetic dipole using 
% VMD or, for a loop source, using empymod.
%
% For empymod, define source center and source dimension of the rectangular
% loop. 
% 
% For VMD, the source is assumed to be located at the origin.
% 
% Specify a regular grid of observation points. The electric field values 
% at each 
% point are provided in EM_field_VMD(n_comp,n_pts,n_t_steps) or 
% EM_field_empymod(n_comp,n_pts,n_t_steps), respectively.

% Create source location and array of observation points distributed 
% regularly on a rectangular Cartesioan grid over a 3D cuboid. 
% Square loop source only defined for empymod.
% Processing sequence: 
% 1. analytic_e_field.m
% 2. plot_single_e_fields.m
% 3. plot_E_fields_threshold_quiver.m
% 4. convolution_vector_3D.m
% 5. plot_single_sensitivities.m
% 6. plot_convolution_threshold_scatter3.m,
% 7. plot_convolution_slice.m,
% 8. plot_convolution_isosurface.m

% Pick analytical solution
tic
Analytical_Solution = pick(2,'VMD','empymod');

switch Analytical_Solution
    case 'empymod'
        fprintf('Analytical solution: %s\n',Analytical_Solution);
        src_center = [0 0 0]; % May be changed for empymod
        src_length2 = pick(1, 0.5, 2.5, 9);
        line = [-src_length2+src_center(1), -src_length2+src_center(2), 0, 1, 1;    % square, centered at [0,0,0]
            src_length2+src_center(1), -src_length2+src_center(2), 0, 1, 1;
            src_length2+src_center(1),  src_length2+src_center(2), 0, 1, 1;
            -src_length2+src_center(1),  src_length2+src_center(2), 0, 1, 1;
            -src_length2+src_center(1), -src_length2+src_center(2), 0, 1, 1];
        src_A = 1/2 * sum(line(1:end-1, 1).*line(2:end, 2) - ...
            line(2:end, 1).*line(1:end-1, 2));
        src_length=2*src_length2;
    case 'VMD'
        src_center = [0 0 0]; % Default for VMD, fixed
        fprintf('Analytical solution: %s\n',Analytical_Solution);
end

% No of components
n_comp = 3;
%n_t_steps=57; % 71 5 dec, 57 4 dec, 43 3 dec

% Set general time vector.
t_exp_min = -7;
t_exp_max = -3;
n_t_steps=(t_exp_max-t_exp_min)*14+1; % 71 5 dec, 57 4 dec, 43 3 dec
t_obs = logspace(t_exp_min, t_exp_max, n_t_steps);

% Conductivity model
sig_air = 1e-9;
sig_earth = 1/10;

% Define grid dimensions and spacings - each node is an observation point
% 
delta_x = 20;
delta_y = 20;
delta_z = 20;
px_min = -40;
px_max = 40;
py_min = -20;
py_max = 20;
pz_min = 0;
pz_max = 20;

% Observation points ( = grid points)
x_coords_vec = px_min:delta_x:px_max;
y_coords_vec = py_min:delta_x:py_max;
z_coords_vec = pz_min:delta_x:pz_max;
% No of points in each dimension
np_x= size(x_coords_vec,2);
np_y= size(y_coords_vec,2);
np_z= size(z_coords_vec,2);
n_pts = np_x*np_y*np_z;
fprintf('Number of observation points: %d %d %d\n',np_x,np_y,np_z);
fprintf('Boundaries (x,y,z): %d %d %d %d %d %d\n',px_min,px_max, ...
    py_min,py_max,pz_min,pz_max);
fprintf('Source length: %d m\n',src_length);


% Format point list
%p_obs = zeros(n_comp,n_pts);
%p_obs[X, Y, Z] = ndgrid(x, y, z);
% Create 3D grid using ndgrid (for volumetric data)
[p_x, p_y, p_z] = ndgrid(x_coords_vec, y_coords_vec, z_coords_vec);

% add 1 mm to points a t the surface to assign them to the conductive earth
% instead of the air space (requirement for empymod)
epsilon = 0%-1e-3;
p_z(:,:,1) = p_z(:,:,1)+epsilon;

% Store points as N×3 list
p_obs = [p_x(:),p_y(:),p_z(:)];

%p_obs = [x_coords_vec,y_coords_vec,z_coords_vec];
%n_p = numel(p_obs);%size(x_coords_vec)*size(y_coords_vec)*size(z_coords_vec);

% Calculate analytical solutions using empymod or VMD
switch Analytical_Solution
    case 'empymod'
        EH = 'E';
        EM_field_empymod = zeros(n_comp,n_pts,n_t_steps);
        % empymod
        % Works but very slow
        E_empymod = app_tem.reference.EH4wire(line(:, 1:3), ...
                p_obs, sig_earth, 0, ...
                'tx_type', 'E', 'rx_type', EH, ...
                't_type', 'time', 't_list', t_obs, ...
                'sigma_air', sig_air);
        EM_field_empymod = E_empymod / src_A; % 3 x n_pts x n_t_steps
        % Store fields and coordinates
        name = ['EM_field_and_coord_empymod_',...
            num2str(np_x),'x',num2str(np_y),'x',num2str(np_z),'_',...
            num2str(t_exp_min),'_',num2str(t_exp_max),...
            '_srcl',num2str(src_length)];
        datafilename = sprintf('%s.mat', name);
        save(datafilename, 'EM_field_empymod', 't_obs', ...
            'x_coords_vec','y_coords_vec','z_coords_vec', ...
            'sig_earth','src_center','src_length');
%{
%empymod
% doesn't work
for ipx = 1:np_x
    for ipy = 1:np_y
        for ipz = 1:np_z
    EH = 'E';
%    title_str3 = ' E';
%    ylabel_str = 'V/Am^3';
    p_obs = [x_coords_vec(ipx),y_coords_vec(ipy),z_coords_vec(ipz)];
    % Get reference solution for finite-length wire src.
    E_empymod = app_tem.reference.EH4wire(line(:, 1:3), ...
                p_obs, sig_earth, 0, ...
                'tx_type', 'E', 'rx_type', EH, ...
                't_type', 'time', 't_list', t_obs, ...
                'sigma_air', sig_air);
    EM_field_ana = E_empymod / src_A;
        end
    end
end
%}

    case 'VMD'
        % VMD
        % Get reference solution for unit dipole src.
            %tmp = @(t) app_tem.examples.homogeneous_hs.getVMDLayeredTransient({'dHz', 'dHr'}, t, ...
            %            norm(point(1:3)), 1/sig_earth, [], 0, 0);
        EM_field_VMD = zeros(n_comp,n_pts,n_t_steps);
        Ez_VMD = deal(zeros(n_t_steps, 1));

        for ip = 1:n_pts
            point(1,1:n_comp) = p_obs(ip,1:n_comp);
            tmp = @(t_obs) app_tem.reference.getVMDLayeredTransient({'Ef'}, t_obs, ...
                norm(point(1,1:2)), 1/sig_earth, [], point(1,3), 0);
            [E_VMD] = deal(zeros(length(t_obs), 1));
            for ii = 1:n_t_steps
                E_VMD(ii) = tmp(t_obs(ii));
            end
            
            EM_field_VMD(1:3,ip,1:n_t_steps) = [-E_VMD*point(1,2)/norm(point(1,1:2)), E_VMD*point(1,1)/norm(point(1,1:2)), Ez_VMD]'; 
            %title_str3 = ' E';
            %ylabel_str = 'V/Am^3';
            %Try this:
            %fprintf('\r\rPoint %6.d out of %6.d',ip, n_pts);
            if ip > 1 fprintf(repmat('\b', 1, 26)); end
            fprintf('Point %6.d out of %6.d',ip, n_pts);
            
        end
        % Store fields and coordinates
        name = ['EM_field_and_coord_VMD_',num2str(np_x),'x',num2str(np_y),...
            'x',num2str(np_z),'_',num2str(t_exp_min),'_',num2str(t_exp_max)];
        datafilename = sprintf('%s.mat', name);
        save(datafilename, 'EM_field_VMD', 't_obs', ...
            'x_coords_vec','y_coords_vec','z_coords_vec', ...
            'sig_earth','src_center');
end
fprintf('\n');
toc
