%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% convolution_vector_3D_split.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab 2024a
%
% Matlab script for performing the 3D convolution of electric fields from 
% a physical and an adjoint source on a regular 3D grid. 
% 
% t_obs and the e-fields have already been created by running
% analytic_e_field.m and written to disk as a mat-file
% (e.g., EM_field_and_coord_empymod_41x21x11.mat).
%
% The physical source is generally a loop source, the adjoint source a VMD
% or a very small loop source. 'split' refers to the fact that both
% physical and adjoint fields are read from separate files.
% Note that one can also use the numerical E-fields instead of the 
% analytical solutions.
%
% The source for the physical forward program is located at the origin. 
% The receiver location for the adjoint approach can be defined. 
%
% Physical field data d1 @ reference point k, transmitter Tx @ i
% Adjoint field data d2 @ reference point k, adjoint source @ Rx receiver 
% location j
%
% Simulated field layout:
%
% S_i,j,k:
%
%      Tx @ i x-----------x Rx @ j
%                   |
%                   |
%                   x
%           reference point @ k
%
% In case Tx and Rx are of the same kind, then for ease of computation 
% observations at the reference point can be
% done for only one Tx. Since the fields are symmetric the coordinates of 
% the reference point are relative to the physical Tx and can be mirrored 
% for the adjoint fields relative to Rx.
%
% Visualization can be done by 
% plot_single_e_and_or_dEdt_fields.m,
% plot_E_fields_threshold_quiver.m,
% plot_single_sensitivities.m, 
% plot_convolution_threshold_scatter3.m,
% plot_convolution_slice.m, and
% plot_convolution_isosurface.m

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define adjoint source center
src_adjoint_center = [20,0,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolve (E and E_adj) or (E and dE_adj/dt)
conv_dEdt = pick(2,false,true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data file for physical fields
% e.g. EM_field_and coord_empymod_41x21x11_-7_-3_srcl5.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Read in data and coordinate mat-file for physical E-fields');
[data_file_phys,location] = uigetfile({'EM_field*.mat'});
if isequal(data_file_phys,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(location,data_file_phys)]);
end

load([location,data_file_phys]);

if contains(data_file_phys, 'VMD')
    disp('VMD physical data read');
    title_string1='VMD';
    EM_field_phys=EM_field_VMD;
    src_length_phys = 0; 
    title_string2_phys=[title_string1,num2str(src_length_phys)];
elseif contains(data_file_phys, 'empymod')
    disp('empymod physical data read');
    title_string1='empymod';
    EM_field_phys=EM_field_empymod;
    src_length_phys = src_length; 
    title_string2_phys=[title_string1,num2str(src_length_phys)];
else
    disp('Input file error physical fields')
    return
end
 

% Detect NaNs in EM_field_phys and set them to 0
EM_field_phys(isnan(EM_field_phys))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data file for adjoint fields
% e.g. EM_field_and coord_empymod_41x21x11_-7_-3.mat
% Note that for a point receiver the data file name should not contain 
% srcl5 or bigger but srcl1 or a VMD
% data_file_phys and data_file_adj may be chosen to be identical 
% for dipole or calculations using small sources 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Read in data and coordinate mat-file for adjoint E-fields');
[data_file_adj,location] = uigetfile({'EM_field*.mat'});
if isequal(data_file_adj,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(location,data_file_adj)]);
end

load([location,data_file_adj]);

if contains(data_file_adj, 'VMD')
    disp('VMD adjoint data read');
    title_string1='VMD';
    EM_field_adj=EM_field_VMD;
    src_length_adj = 0; 
    title_string2_adj=[title_string1,num2str(src_length_adj)];
elseif contains(data_file_phys, 'empymod')
    disp('empymod adjoint data read');
    title_string1='empymod';
    EM_field_adj=EM_field_empymod;
    src_length_adj = src_length; 
    title_string2_adj=[title_string1,num2str(src_length_adj)];
else
    disp('Input file error adjoint fields')
    return

end

% Detect NaNs in EM_field and set them to 0
EM_field_adj(isnan(EM_field_adj))=0;

% Convolve (E and E_adj) or (E and dE_adj/dt)
if conv_dEdt
    fprintf('Convolution over dEdt!\n')
else
    fprintf('Convolution over E!\n')
end

% Sanity check
if ~isequal(size(EM_field_phys),size(EM_field_adj))
    disp('Dimensions of physical and adjoint fields do not agree!');
    return
end

% Original size parameters
Nx = numel(x_coords_vec);  % 41
Ny = numel(y_coords_vec);  % 21
Nz = numel(z_coords_vec);  % 11
n_pts = Nx * Ny * Nz;      % 9471
n_comp = size(EM_field_phys, 1);   % 3
n_t_obs = size(EM_field_phys, 3);  % 43

% Sanity check
if ~isequal(n_pts-size(EM_field_phys, 2),0)
    disp('Dimensions of points and fields do not agree!');
    return
end
fprintf('Source location: %d %d %d\n',src_center);
fprintf('Source length physical source: %d m\n',src_length_phys);
fprintf('Source length adjoint source: %d m\n',src_length_adj);

% Generate original grid (vectorized 3D locations)
[X, Y, Z] = ndgrid(x_coords_vec, y_coords_vec, z_coords_vec);  % size: [Nx, Ny, Nz]
X = X(:); Y = Y(:); Z = Z(:);  % [9471 x 1] column vectors

% Create map from (x, y, z) to linear index
% Use a containers.Map for fast lookup
key_from_coords = @(x, y, z) sprintf('%.1f_%.1f_%.1f', x, y, z);
coord_to_index = containers.Map();
for k = 1:n_pts
    key = key_from_coords(X(k), Y(k), Z(k));
    coord_to_index(key) = k;
end

% Find all valid receiver indices and their matching sources
valid_receiver_indices = [];
source_indices = [];

% Define Rx location for adjoint field source
x_adjoint_src = src_adjoint_center(1);
y_adjoint_src = src_adjoint_center(2);
z_adjoint_src = src_adjoint_center(3);
fprintf('Receiver location: %d %d %d\n',src_adjoint_center);

% Loop through all points that have a valid left-shifted partner
for k = 1:n_pts
    x = X(k); y = Y(k); z = Z(k);
    x_shifted = x - x_adjoint_src;
    y_shifted = y - y_adjoint_src;
    z_shifted = z - z_adjoint_src;

    % Check if the shifted point is still inside the domain
    key_shifted = key_from_coords(x_shifted, y_shifted, z_shifted);
    if isKey(coord_to_index, key_shifted)
        valid_receiver_indices(end+1) = k;
        source_indices(end+1) = coord_to_index(key_shifted);
    end
end
fprintf('Total number of grid points: %d\n',n_pts);
fprintf('Valid reference points: %d\n',size(valid_receiver_indices,2));

n_valid_pts = numel(valid_receiver_indices);
n_uniform = 40000;
conv_nonuniform_matrix = zeros(n_t_obs, n_valid_pts);  % Each column is a time series
%conv_linear_matrix = zeros(2*n_uniform - 1, n_valid_pts);
conv_spline_matrix = zeros(2*n_uniform - 1, n_valid_pts);

t_uniform = linspace(t_obs(1), t_obs(end), n_uniform); % fine uniform time grid
t_min = min(t_obs);
t_max = max(t_obs);
    
% Loop over valid receiver-source pairs

for k = 1:n_valid_pts
    ip1 = source_indices(k);            % source index (shifted x)
    ip2 = valid_receiver_indices(k);    % receiver index

    d1 = squeeze(EM_field_phys(:,ip1,:));
    d2 = squeeze(EM_field_adj(:,ip2,:));
    %n_comp = size(d1,1); % No of components, usually 3.

    % --- Step 1: Non-uniform convolution ---
    % Non-uniform convolution
    conv_nonuniform = zeros(1, n_t_obs);
    f_tau = zeros(n_comp,n_t_obs);

    % Loop through non-uniform time steps
    if k == 1 
        fprintf('Execute nonuniform convolution\n'); 
    end
    for it = 1:n_t_obs
        tau = t_obs(1:it);  % Time points up to current
        if conv_dEdt
            d1_temp = d1(:,1:it);
            f_tau = gradient(d1_temp,tau,2);   % Corresponding values of d1
        else
            f_tau = d1(:,1:it);   % Corresponding values of d1
        end
        g_t_minus_tau = zeros(n_comp, it);

        % g(t - tau): Interpolation of d2 at shifted time points t(i) - tau
        for ic = 1:n_comp
            %g_interp = interp1(t_obs, d2(ic, :), t_obs(it) - tau, 'linear', 'extrap');
            g_interp = interp1(t_obs, d2(ic, :), t_obs(it) - tau, 'linear', 0);
            g_t_minus_tau(ic,:) = g_interp; 
        end

        integrand = sum(f_tau .* g_t_minus_tau,1);  % Pointwise multiplication
    
        if numel(tau) > 1  % Only integrate if there is more than one point
            conv_nonuniform(it) = trapz(tau, integrand);
        else
            conv_nonuniform(it) = 0;  % No integration needed for a single point
        end
    end
    conv_nonuniform_matrix(:,k) = conv_nonuniform.';  

    % --- Step 3: Uniform grid (spline resampling) ---
    dt = t_uniform(2)-t_uniform(1); % assume equidistant time points
    
    if k == 1 
        fprintf('Execute spline convolution %d time steps\n',n_uniform);
        fprintf('Point %6.d out of %6.d',k, n_valid_pts);
    end
    d1_spline = zeros(n_comp,n_uniform);
    d2_spline = zeros(n_comp,n_uniform);

    for ic = 1:n_comp
        d_temp = interp1(t_obs, d2(ic,:), t_uniform, 'spline');
        d2_spline(ic,:) = d_temp;

        if conv_dEdt
            d_interp = interp1(t_obs, d1(ic,:), t_uniform, 'spline');
            d_temp = gradient(d_interp,t_uniform,2);
        else
            d_temp = interp1(t_obs, d1(ic,:), t_uniform, 'spline');
        end
        d1_spline(ic,:) = d_temp;
    end

    conv_spline = zeros(1, 2*n_uniform - 1);

    for ic = 1:n_comp
        c = conv(d1_spline(ic,:), d2_spline(ic,:)) * dt; % scale by dt
        conv_spline = conv_spline + c;
        %conv_linear = conv_linear(1:length(t_uniform)); % same length
    end
    conv_spline_matrix(:,k) = conv_spline;
    if k > 1 
        fprintf(repmat('\b', 1, 20)); 
        fprintf('%6.d out of %6.d',k, n_valid_pts);
    end
end % end loop over valid receiver locations

% Store convolution results and coordinates
if conv_dEdt
        name = ['Convolution_Results_',title_string2_phys,'_',title_string2_adj,...
            '_grad_',num2str(Nx),'x',num2str(Ny),...
            'x',num2str(Nz),'_',num2str(x_adjoint_src),'_',...
            num2str(y_adjoint_src),'_',num2str(z_adjoint_src),'_',...
            num2str(log10(t_min)),'_',num2str(log10(t_max))];
else
        name = ['Convolution_Results_',title_string2_phys,'_',title_string2_adj,...
            '_',num2str(Nx),'x',num2str(Ny),...
            'x',num2str(Nz),'_',num2str(x_adjoint_src),'_',...
            num2str(y_adjoint_src),'_',num2str(z_adjoint_src),'_',...
            num2str(log10(t_min)),'_',num2str(log10(t_max))];
end
        
datafilename = sprintf('%s.mat', name);
fprintf('\nSaving data file: %s\n',datafilename);

if strcmp(title_string2_phys,'VMD0') % VMD0_phys has always VMD0_adj
        EM_field_VMD_phys = EM_field_phys;
        EM_field_VMD_adj  = EM_field_adj;
        save(datafilename, 'conv_nonuniform_matrix', 't_obs', ...
            'conv_spline_matrix','t_uniform',...
            'x_coords_vec','y_coords_vec','z_coords_vec', ...
            'sig_earth', ...
            'src_center','src_length_phys', ...
            'src_adjoint_center','src_length_adj', ...
            'EM_field_VMD_phys','EM_field_VMD_adj', ...
            'valid_receiver_indices','-v7.3');
elseif (contains(title_string2_phys,'empymod') && contains(title_string2_adj,'VMD'))
        EM_field_empymod_phys = EM_field_phys;
        EM_field_VMD_adj  = EM_field_adj;
        save(datafilename, 'conv_nonuniform_matrix', 't_obs', ...
            'conv_spline_matrix','t_uniform',...
            'x_coords_vec','y_coords_vec','z_coords_vec', ...
            'sig_earth', ...
            'src_center','src_length_phys', ...
            'src_adjoint_center','src_length_adj', ...
            'EM_field_empymod_phys','EM_field_VMD_adj', ...
            'valid_receiver_indices','-v7.3');
elseif (contains(title_string2_phys,'empymod') && contains(title_string2_adj,'empymod'))
        EM_field_empymod_phys = EM_field_phys;
        EM_field_empymod_adj  = EM_field_adj;
        save(datafilename, 'conv_nonuniform_matrix', 't_obs', ...
            'conv_spline_matrix','t_uniform',...
            'x_coords_vec','y_coords_vec','z_coords_vec', ...
            'sig_earth', ...
            'src_center','src_length_phys', ...
            'src_adjoint_center','src_length_adj', ...
            'EM_field_empymod_phys','EM_field_empymod_adj', ...
            'valid_receiver_indices','-v7.3');

end
fprintf('%s written\nDone\n ',datafilename);
        
toc