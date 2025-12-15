function F = EH4wire(tx, rx, sigma, depth, varargin)
    % (Semi)analytic E-/H-field for finite bipole at surface of layered HS.
    %
    % The routine is a wrapper to empymod by D.Werthmüller using
    % Matlab-Python-API
    % https://github.com/emsig/empymod
    %
    % SYNTAX
    %   F = EH4wire(tx, rx, sigma, depth, [varargin])
    %
    % INPUT PARAMETER
    %   tx ... Matrix [n x 3] of wire start and end points
    %          Note: last end point will be handled as next
    %                starting point, such that for close loops the initial
    %                point should occur twice!
    %   rx ... Matrix [n x 3] of observation points
    %          Note: tx_z and rx_z have to be equal for all points
    %   depth   ... Vector of layer boundaries (length(sigma))
    %   sigma   ... Vector of layer conductivities
    %
    % OPTIONAL PARAMETER
    %   t_type  ... Character denoting if default: 'time' oder 'freq'uency
    %               domain problem is considered
    %               Note: if time -> impulse time-domain response
    %   t_list  ... Vector of times or frequencies ['time'; 'freq']
    %   tx_type ... Character denoting if 'H'magnetic or
    %               default: 'E'electric source is considered
    %   rx_type ... Character [default: 'dH', 'H', 'dE', 'E] denoting if 
    %               magnetic or electric field (time derivative) should be 
    %               calculated
    %   epsilon ... Vector of relative horz./vert. magnetic permeabilities
    %               Note: we set  -> horz. = vert.
    %   sigma_air ... Scalar: air conductivity [default 1e-14]
    %   srcpts  ... int: Number of integration points for bipole tx
    %                         < 3  : bipole, calculated as dipole at centre
    %                         >= 3 : bipole
    %   I       ... Scalar: Source current strength
    %                       0, output normalized to tx of 1m length,
    %                          and tx strength of 1A
    %                    != 0, output is returned for given tx length,
    %                          and source strength
    %   verb    ... Scalar, level of verbosity {0, 1, 2, 3, 4}
    %
    % OUTPUT PARAMETER
    %   F ... Matrix [3 x n_rx x n_t] of fields at rx positions
    %
    % REQUIREMENTS / REMARKS
    %   Tested with Matlab 2024b, Pyhton 3.10.12, numpy 1.22.4,
    %               scipy 1.13.1, empymod 2.2.0
    %   By default: solution for TEM problem is considered.

    %% Check and fetch input.

    assert(ismatrix(tx) && size(tx, 2) == 3 && size(tx, 1) >= 2);
    assert(ismatrix(rx) && size(rx, 2) == 3);
    assert(isvector(sigma) && length(sigma) >= 1);
    assert(isvector(depth) && length(depth) == length(sigma));

     % Define possible input keys and its properties checks.
    input_keys = {'tx_type', 'rx_type', 't_type', 't_list', ...
                  'epsilon', 'I', 'srcpts', 'verb', 'sigma_air'};

    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, 'E', @(x) ischar(x) && isscalar(x));
    parser_obj.addParameter(input_keys{2}, 'dH', @(x) ischar(x) && (isscalar(x) || length(x) == 2));
    parser_obj.addParameter(input_keys{3}, 'time', @ischar);
    parser_obj.addParameter(input_keys{4}, 1, @isvector);
    parser_obj.addParameter(input_keys{5}, [], @(x) isvector(x) && length(x) == length(sigma));
    parser_obj.addParameter(input_keys{6}, 1, @isscalar);
    parser_obj.addParameter(input_keys{7}, 12, @(x) isscalar(x) && x > 0);
    parser_obj.addParameter(input_keys{8}, 1, @(x) isscalar(x) && x < 5 && x >= 0);
    parser_obj.addParameter(input_keys{9}, 1e-14, @(x) isscalar(x) && x > 0);

    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    % Transform input parameter.
    % Physical parameter.
    rho = 1./[args.sigma_air, sigma(:).'];
    if ~isempty(args.epsilon)
        epsilon = [0, args.epsilon(:).'];
    else
        epsilon = [0, ones(size(sigma(:).'))];
    end

    % Source domain type.
    n_t = length(args.t_list);
    switch args.t_type
        % Time domain.
        case 'time'
            % Define solution type.
            if isscalar(args.rx_type)
                % 'Switch-off response' i.e. field itselfe
                signal = -1;
            else
                % 'Impulse response' i.e.    field time-derivative
                signal = 0;
            end
        % Frequency domain.
        % -> Solution is given in f-domain too.
        case 'freq'
            signal = py.None;
        otherwise
            error('Unsupported domain type.');
    end

    % RX/TX type
    switch args.tx_type
        case 'H'
            msrc = true;
        case 'E'
            msrc = false;
        otherwise
            error('Unsupported TX type.');
    end
    switch args.rx_type(end)
        case 'H'
            mrec = true;
        case 'E'
            mrec = false;
        otherwise
            error('Unsupported RX type.');
    end

    % Set source as list of bipoles: [x0, x1, y0, y1, z0, z1]
    n_tx_wire = size(tx, 1)-1;
    py_tx = zeros(n_tx_wire, 6);
    for ii = 1:n_tx_wire
        py_tx(ii, :) = reshape([tx(ii,:);tx(ii+1,:)], [], 6);
    end

    % Receiver points.
    n_rx = size(rx, 1);
    % Set rx as array of dipoles: [{x}, {y}, {z}].
    py_rx = {rx(:, 1).', rx(:, 2).', rx(:, 3).'};

	%% Do the calculation.

    % Set rx azimuth, dip for R3.
    % azimuth (°) ... horizontal deviation from x-axis, anti-clockwise.
    % dip (°)     ... vertical deviation from xy-plane downwards.
    azi_dip = {0,  0;  % x
               90, 0;  % y
               0, 90}; % z

    % Collect arguments.
    args_fix = {'epermH', epsilon, 'epermV', epsilon, 'verb', args.verb, ...
                'freqtime', args.t_list, 'signal', signal, 'depth', depth, ...
                'res', rho, 'srcpts', args.srcpts, 'strength', args.I, ...
                'mrec', mrec, 'msrc', msrc};

    % Loop over wire segments.
    [F, F_] = deal(zeros(3, n_rx, n_t));
    for jj = 1:n_tx_wire
        fprintf('wire %d/%d\n',jj,n_tx_wire);
        args_tx = {'src', py_tx(jj, :)};
        % Loop over receiver orientations.
        for ii = 1:size(azi_dip, 1)
            fprintf('azi %d/%d\n',ii,size(azi_dip, 1))
            args_rx = {'rec', [py_rx, repmat(azi_dip{ii, 1}, 1, n_rx), ...
                                      repmat(azi_dip{ii, 2}, 1, n_rx)]};
            F_(ii,:,:) = py.empymod.bipole(pyargs(args_fix{:}, args_tx{:}, args_rx{:})).';
        end
        % Sum up separate wire solutions.
        F = F + F_;
    end
end
