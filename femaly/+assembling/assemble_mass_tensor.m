function T = assemble_mass_tensor(dofmap)
    % Assemble operator \int \sigma_m \phi_i \phi_j
    %
    %   \phi_i       - H^1-conforming basis functions
    %   \sigma_m     - P_0 conforming basis functions
    %
    % INPUT PARAMETER
    %   dofmap ... Struct, containing mesh and FE element objects
    %              (the latter contain H^1 conforming basis functions)
    %              as well as the cell-2-DOF mapping.
    %
    % OUTPUT PARAMETER
    %   T ... Sparse 3-way-tensor [dofmap.dim x dofmap.dim], representing
    %         mass matrix fderivative w.r.t. P_0 parameter.

    local_element_dim = dofmap.element.fe_space_dim;
    dim = dofmap.mesh.dim;
    num_cells = size(dofmap.mesh.cells, 2);

    % Pick quadrature rule
    quad_degree = 2*dofmap.element.order;
    [x, w] = fe.make_quadrature_rule(dim, quad_degree);
    num_quad_points = size(x, 1);

    % Tabulate basis and basis curls at quadrature points
    basis_curls = zeros(local_element_dim, dofmap.element.curl_shape, num_quad_points);
    basis       = zeros(local_element_dim, dim, num_quad_points);
    for k = 1:num_quad_points
        basis_curls(:, :, k) = dofmap.element.tabulate_basis_curl(x(k, :));
        basis      (:, :, k) = dofmap.element.tabulate_basis     (x(k, :));
    end

    % Fetch some data and preallocate temporaries
    cells = dofmap.mesh.cells;
    coords = dofmap.mesh.vertex_coords;
    cell_dofs = dofmap.cell_dofs;
    A = zeros(local_element_dim, local_element_dim);
    jac = zeros(dim, dim);
    jac_inv = zeros(dim, dim);
    temp = zeros(local_element_dim, dim);
    detJ = zeros(1, 1, 'double');
    detJinv = zeros(1, 1, 'double');

    % Preallocate assembly data
    nnz = num_cells*local_element_dim^2;
    I = zeros(nnz, 1, 'double');
    J = zeros(nnz, 1, 'double');
    V = zeros(nnz, 1, 'double');
    offsets = uint32(1:local_element_dim^2);

    % Loop over cells
    for c = 1:num_cells

        % Compute geometric quantities
        jac(:, :) = coords(:, cells(1:dim, c)) - coords(:, cells(dim+1, c));
        detJ(1) = abs(det(jac));
        detJinv(1) = 1/detJ;
        jac_inv(:, :) = inv(jac);

        % Zero cell matrix
        A(:, :) = 0;

        % Loop over quadrature points and assemble cell integrals
        for k = 1:num_quad_points
            temp(:, :) = basis(:, :, k)*jac_inv;
            A(:, :) = A + detJ*w(k)*(temp*temp.');
        end

        % Compute global dof indices and store cell matrix
        [J(offsets), I(offsets)] = meshgrid(cell_dofs(:, c));
        V(offsets) = A;
        offsets(:) = offsets + local_element_dim^2;

    end

    % Create sparse 3-way-tensor from index and value vectors.
    size_M = [dofmap.dim, dofmap.dim, num_cells];
    index_M = [I, J, kron((1:num_cells).', ones(local_element_dim^2, 1))];
    T = sptensor.Tensor3Coord(size_M, index_M, V);
end
