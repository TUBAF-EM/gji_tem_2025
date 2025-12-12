function sol = assemble(type, mesh, FE_order, bnd_sum, tem_info)
    % Provides handles to matrices, rhs and FE info for TEM.
    %
    % SYNTAX
    %   sol = assemble(type, mesh, FE_order, bnd_sum, tem_info)
    %
    % INPUT PARAMETER
    %   type     ... Character, denoting the problem dimension [3D].
    %   mesh     ... Object of mesh class.
    %   FE_order ... Scalar, denoting finite element oder.
    %   bnd_sum  ... Cell [n x 2], containing:
    %               bnd_sum{n, 1}
    %                   character denoting the BC type
    %               bnd_sum{n, 2} containing either:
    %                   Cell [m x k = 1...2]: equal sized cells containing
    %                       k = 1: function handle(s) to identify where bc
    %                              (parts) should be applied
    %                       k = 2: function handle(s) taking coordinate and
    %                              providing function value to be applied
    %                              at bc (parts)
    %               or
    %                   Cell [m x k = 1...3]:
    %                       k = 1: Vector of facet tags/ids of every facet
    %                              in mesh
    %                       k = 2: Cell(s) containing vector of all
    %                              boundary facets belonging to considered
    %                              bc (parts)
    %                       k = 2: Cell(s) containing function handle
    %                              taking coordinate and providing function
    %                              value to be applied at bc (parts)
    %   tem_info ... Struct, containing survey info and mappings to obtain
    %                observation from simulated quantities.
    %
    % OUTPUT PARAMETER
    %   sol        ... Struct, containing:
    %   O          ... Matrix (observation/measurement operator).
    %   dofmap     ... Struct, containing cell2DOF mapping as well as the
    %                  underlying mesh and FE element objects.
    %   type       ... Character, denoting problem dimension / type.
    %   KMb_handle ... Handle taking vector 'param' providing
    %                   K - curl-curl matrix
    %                   M - mass matrix (depending on model parameter)
    %                   b - Vector (rhs)
    %   TM         ... Mass matrix tensor representation as sptensor object

    % Fetch and check input.
    assert(ischar(type) && any(str2double(type(1)) == 3));
    assert(isa(mesh, 'meshing.Mesh'));
    assert(isscalar(FE_order));
    assert(iscell(bnd_sum) && size(bnd_sum, 2) == 2);
    assert(strcmp(bnd_sum{1}, 'dirichlet'));
    assert(isstruct(tem_info) && all(isfield(tem_info, {'pt_TX', 'lm', 'srcm'})));
    if isfield(tem_info, 'data_type')
        assert(iscell(tem_info.data_type) && ...
               all(cellfun(@(x) any(strcmp({'dHdt', 'dBdt', 'E'}, x)), ...
                       tem_info.data_type)));
        assert(length(unique(tem_info.data_type)) == length(tem_info.data_type)); % avoid multiples
    else
        tem_info.data_type = {'dBdt'};
    end
    %
    sol = struct();
    sol.type = [num2str(mesh.dim),'D'];

    % Build Nedelec function space.
    element = fe.create_nedelec_element(mesh.dim, FE_order);
    dofmap = assembling.build_dofmap(mesh, element);
    sol.dofmap = dofmap;

    % Loop over source parts.
    assert((isscalar(tem_info.srcm) && size(tem_info.pt_TX, 2) == 3) || ...
           (length(tem_info.srcm) == length(unique(tem_info.pt_TX(:, end)))) );
    tmp_idx = max(tem_info.srcm);
    tmp_pt_m = unique(tem_info.pt_TX(:, end));
    b = zeros(dofmap.dim, 1);
    for ll = 1:length(tem_info.srcm)
        cur_pt_TX = tem_info.pt_TX(tem_info.pt_TX(:, end) == tmp_pt_m(ll), :);
        if isequal(cur_pt_TX(1, 1:3), cur_pt_TX(end, 1:3))
            % Loop area from 2D polygone corner points (wiki).
            scale = 1/2 * sum(cur_pt_TX(1:end-1, 1) .* cur_pt_TX(2:end, 2) - ...
                              cur_pt_TX(2:end, 1)   .* cur_pt_TX(1:end-1, 2));
            if ~all(cur_pt_TX(1, 3) == cur_pt_TX(2:end, 3))
                warning(['Source area is derived from projection of points on xy-plane.', ...
                         ' Deviations to true area due to topography expected.']);
            end
        else
            % Length of wire from points.
            tmp_length = 0;
            for ww = 1:size(cur_pt_TX, 1)-1
                tmp_length = tmp_length + vecnorm(cur_pt_TX(ww+1, 1:3) - cur_pt_TX(ww, 1:3));
            end
            scale = tmp_length;
        end

        % Ensure current wire mesh edges to have same orientation.
        % -> Partition (tx) wire into plus- and minus-oriented edges.
        assert(any(unique(tem_info.lm) == tem_info.srcm(ll)), ...
               'No line elements from desired markers found in mesh.');
        oriented_lm = meshing.orient_wire(mesh, tem_info.lm, tem_info.srcm(ll), ...
                                          10*tmp_idx+ll, 20*tmp_idx+ll);
        % Assemble right-hand side, scaled to unit strength.
        b = b + assembling.assemble_edge_source(dofmap, oriented_lm, 10*tmp_idx+ll, ...
            @(x, tau) 1/scale, 0);
        b = b - assembling.assemble_edge_source(dofmap, oriented_lm, 20*tmp_idx+ll, ...
            @(x, tau) 1/scale, 0);
    end

    % Homogeneous Dirichlet on boundary everywhere.
    [bdofs, ~] = assembling.build_dirichlet_dofs(dofmap, bnd_sum{2}{:});

    % Assemble sigma-independing operators.
    % Curl-Curl matrix.
    K = assembling.assemble_curl_curl(dofmap, 0, util.EMConstants.mu_0);
    % Apply DBC.
    [K, b] = assembling.apply_dirichlet_bc(K, b, bdofs, 0);
    % Mass matrix derivative tensor.
    TM = assembling.assemble_mass_tensor(dofmap);
    sol.TM = TM; % required for sensitivity calculation

    if isfield(tem_info, 'pt_RX')
        % Assembling measurement operator.
        O = [];
        for cc = 1:length(tem_info.data_type)
            switch tem_info.data_type{cc}
                case {'dBdt', 'dHdt'}
                O_ = assembling.assemble_point_sources(dofmap, tem_info.pt_RX(:,1:3).', 'curl');
                if strcmp(tem_info.data_type, 'dHdt')
                    O_ = O_ ./ util.EMConstants.mu_0;
                end

                case {'E'}
                O_ = assembling.assemble_point_sources(dofmap, tem_info.pt_RX(:,1:3).');
            end

            % Vertically concatenate different data types.
            O = [O, O_]; %#ok<AGROW> only few different types exected
        end
    else
        % d = u, i.e. the FE coefficient vector.
        O = speye(sol.dofmap.dim);
    end
    sol.O = O;
    sol.KMb_handle = @(sigma) assemble_system(sigma, K, b);

    % Helper.
    function [K, M, b] = assemble_system(sigma, K, b)

        % Assemble sigma-dependent operators.
        M_ = ttv(TM, sigma, 3);
        M = assembling.apply_dirichlet_bc(M_, [], bdofs, 0);
    end
end
