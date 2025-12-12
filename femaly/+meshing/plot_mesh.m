function plot_mesh(mesh, varargin)
    % Plot mesh, optionally with cell, facet, edge, and/or vertex markers
    %
    % SYNTAX
    %   plot_mesh(mesh[, 'Name', Value, ...])
    %
    % INPUT PARAMETER
    %   mesh           ... Object of Mesh class
    %
    % OPTIONAL PARAMETERS (NAME-VALUE PAIRS)
    %   cell_markers   ... Vector [num_cells x 1]
    %   facet_markers  ... (sparse) vector [num_facets x 1]
    %   edge_markers   ... (sparse) vector [num_edges x 1]
    %   vertex_markers ... (sparse) vector [num_vertices x 1]
    %   plane          ... visualizes a slice through the cells
    %                      intersecting a plane given by the implicit
    %                      equation:
    %
    %                          plane(x(1:3, :)) == 0;
    %
    %                      'plane' is a function handle of signature:
    %
    %                          double = function(double(1:3, :));
    %
    %                      for example, xz-plane intersecting point [0; 1/2; 0]:
    %
    %                          plane = @(x) 0*x(1,:) + 1*x(2,:) + 0*x(3,:) - 1/2;
    %
    %                      the function shall be vectorized along the second
    %                      axis of x.
    %                      If omitted, xz-plane intersecting point [0; 0; 0]
    %                      is considered.
    %
    % NOTES
    %   Facet, edge, and vertex markers will be typically sparse
    %   vectors with integer values but this is not enforced.
    %   In particular, zeros are ignored rather than plotted
    %   and legend features labels with all (integer) values.
    %
    %   On the other hand, cell markers are interpreted as
    %   continuous, double-valued field and it is accompanied
    %   by a color bar rather than a legend with values.
    %
    % EXAMPLE (2D)
    %   mesh = meshing.generate_mesh2D();
    %   mesh.compute_connectivity(1, 0);
    %   cm = mod(1:mesh.num_entities(mesh.dim), 3);
    %   em = mod(1:mesh.num_entities(1), 4);
    %   vm = mod(1:mesh.num_entities(0), 5);
    %   meshing.plot_mesh(mesh, 'cell_markers', cm, 'edge_markers', em, ...
    %                     'vertex_markers', vm);
    %
    % EXAMPLE (3D)
    %   mesh = meshing.generate_mesh3D();
    %   mesh.compute_connectivity(1, 0);
    %   cm = mod(1:mesh.num_entities(mesh.dim), 3);
    %   em = mod(1:mesh.num_entities(1), 4);
    %   vm = mod(1:mesh.num_entities(0), 5);
    %   pl = @(x) 3*x(1,:) + 1*x(2,:) + 2*x(3,:) - 2;
    %   meshing.plot_mesh(mesh, 'cell_markers', cm, 'edge_markers', em, ...
    %                     'vertex_markers', vm, 'plane', pl);
    %
    % EXAMPLE (3D)
    %   mesh = meshing.generate_unit_cube_mesh([3, 3, 3]);
    %   pl = @(x) x(2,:) - 1/3;
    %   meshing.plot_mesh(mesh, 'plane', pl, 'preferred_side', +1);

    assert(isa(mesh, 'meshing.Mesh'));

    isvec = @(x) isempty(x) || isvector(x) && isnumeric(x) && isreal(x);
    isfun = @(x) isempty(x) || isa(x, 'function_handle');

    parser_obj = inputParser();
    parser_obj.addParameter('cell_markers',   [], isvec);
    parser_obj.addParameter('facet_markers',  [], isvec);
    parser_obj.addParameter('edge_markers',   [], isvec);
    parser_obj.addParameter('vertex_markers', [], isvec);
    parser_obj.addParameter('plane',          [], isfun);
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    if ~isempty(args.facet_markers) && mesh.dim == 3
        error('plotting of tetrahedra facet markers not yet implemented!');
    end

    % Facets are edges in 2D
    if ~isempty(args.facet_markers) && mesh.dim == 2
        if ~isempty(args.edge_markers)
            error('supply only one of "facet_markers" and "edge_markers" arguments!');
        else
            args.edge_markers = args.facet_markers;
        end
    end

    if ~isempty(args.plane) && mesh.dim ~= 3
        error('intersection with plane only implemented in 3D!');
    elseif isempty(args.plane) && mesh.dim == 3
        args.plane = @(x) 0*x(1,:) + 1*x(2,:) + 0*x(3,:);
    end

    hold on;
    axis equal;
    axis tight;

    annotated_plots = [];

    % Plot mesh and cell markers (if given)
    switch mesh.dim
    case 2
        assert(isempty(args.plane));
        plot_cell_markers_2d(mesh, args.cell_markers);
    case 3
        plot_cell_markers_3d(mesh, args.cell_markers, args.plane);
    end

    % Plot edge markers if given
    if ~isempty(args.edge_markers)
        plots = plot_edge_markers(mesh, args.edge_markers);
        annotated_plots = [annotated_plots(:); plots(:)];
    end

    % Plot vertex markers if given
    if ~isempty(args.vertex_markers)
        plots = plot_vertex_markers(mesh, args.vertex_markers);
        annotated_plots = [annotated_plots(:); plots(:)];
    end

    % Issue legend if appropriate
    if ~isempty(annotated_plots)
        legend(annotated_plots, 'Location', 'BestOutside');
    end

    hold off;
end


function plot_cell_markers_2d(mesh, cell_markers)
    if isempty(cell_markers)
        coords = mesh.vertex_coords;
        triplot(double(mesh.cells.'), coords(1, :), ...
                coords(2, :), 'Color', 'k');
    else
        patch('Faces', mesh.cells.', ...
              'Vertices', mesh.vertex_coords.', ...
              'FaceVertexCData', cell_markers(:), ...
              'FaceColor', 'flat', ...
              'EdgeColor', 'none');
        colorbar();
    end
end


function plot_cell_markers_3d(mesh, cell_markers, plane)

    % Define test if cell contains vertices on both side of the plane.
    intersects_ = @(sdist) any(sdist>=0) && any(sdist<=0);

    % Compute the normal vector of the plane.
    vec_normal = plane(eye(3)) - plane(zeros(3));
    % Find two orthogonal vectors to define the local 2D plane.
    [~, ~, V] = svd(vec_normal(:)');
    vec_orth_1 = V(:, 2); % First orthogonal vector
    vec_orth_2 = V(:, 3); % Second orthogonal vector

    % Loop over all cells.
    num_cells = mesh.num_entities(mesh.dim);
    intersected_cells = false(1, num_cells);
    max_intersection_cells = size(mesh.cells, 2);                         % Each tetrahedron can contribute up to 4 intersection points
    intersection_points = zeros(max_intersection_cells * 4, 3);           % Preallocate for 3D points
    % intersection_points_projected = zeros(max_intersection_cells * 4, 2); % Preallocate for 2D points
    intersection_faces = NaN(max_intersection_cells, 4);                  % Preallocate for faces
    point_count = 0;
    face_count = 0;
    for c = 1:num_cells

        % Get current cell vertex coords and signed distance of each vertex
        % in cell to the plane.
        cell_coords = mesh.vertex_coords(:, mesh.cells(:, c).');
        tmp = plane(cell_coords);

        % Check, if cell intersects plane if there are vertices on both
        % sides.
        if ~intersects_(tmp)
            continue
        end

        % Find the intersection points of the cell edges with the plane.
        points = zeros(4, 3); % Preallocate for up to 4 intersection points per tetrahedron
        point_idx = 0;        % Local counter for points in this tetrahedron
        for j = 1:4
            for k = j+1:4
                % Get the two vertices of the edge and their distance to
                % the plane.
                v1 = cell_coords(:, j);
                v2 = cell_coords(:, k);
                tmp = plane([v1, v2]);

                % Check if the edge crosses the slicing plane.
                if intersects_(tmp)
                    % Compute the parametric line coordinate t.
                    if all(tmp == 0)
                        t = 0;
                    else
                        t = tmp(1) / (tmp(1) - tmp(2));
                    end

                    % Compute the intersection point.
                    intersection_point = (1 - t) * v1 + t * v2;
                    point_idx = point_idx + 1;
                    points(point_idx, :) = intersection_point;
                end
            end
        end

        % Trim unused/multiple rows in the local points array.
        points = points(1:point_idx, :);
        points = unique(points, 'rows', 'stable');
        if isempty(points) || size(points, 1) < 3
            continue
        end

        % Project points onto the local 2D plane.
        projected_points = [points * vec_orth_1, points * vec_orth_2];

        % If there are at least 3 intersection points, form a polygon.
        if size(points, 1) == 3
            % Default triangle ordering.
            tri = [1, 2, 3];

        elseif size(points, 1) == 4

        % Find centroid and angles.
        center = mean(projected_points, 1);
        a = atan2(projected_points(:, 2) - center(2), projected_points(:, 1) - center(1));

        % Sort quadrilateral points clockwise.
        [~, tri] = sort(a);

        else
            error('Expeting less than 4 points.');
        end

        % Store the intersection points and faces.
        num_new_points = size(points, 1);
        % intersection_points_projected(point_count+1:point_count+num_new_points, :) = projected_points;
        intersection_points(point_count+1:point_count+num_new_points, :) = points;
        point_count = point_count + num_new_points;

        % Store polygons.
        tmp = point_count-length(tri)+1:point_count;
        intersection_faces(face_count+1, 1:length(tri)) = tmp(tri);
        face_count = face_count + 1;
        intersected_cells(c) = true; % Store cell idx
    end

    % Trim unused rows in the preallocated arrays.
% TODO: - Speed up by restricting to unique points?
%       - Fix: intersection_points are not 100% coplanar with the plane
%              such that patch plot is actually not 2D
    % intersection_points_projected = intersection_points_projected(1:point_count, :);
    intersection_points = intersection_points(1:point_count, :);
    intersection_faces = intersection_faces(1:face_count, :);

    % Plot the 2D intersection surfaces in 3D and rotated appropriately.
    if isempty(cell_markers)
        patch('Faces', intersection_faces, 'Vertices', intersection_points, ...
              'FaceColor', 'cyan', 'EdgeColor', 'black', 'EdgeAlpha', 0.25);
    else
        cell_markers = cell_markers(:);
        patch('Faces', intersection_faces, 'Vertices', intersection_points, ...
              'FaceVertexCData', cell_markers(intersected_cells), ...
              'FaceColor', 'flat', 'EdgeColor', 'black', 'EdgeAlpha', 0.25);
        colorbar;
    end
    view(-vec_normal);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
end


function plts = plot_edge_markers(mesh, edge_markers)

    coords = mesh.vertex_coords;
    edges_to_vertices = mesh.get_connectivity(1, 0);

    % Build unique marker values and mapping back to marked edges
    [marked_edges, ~, values] = find(edge_markers(:));
    [values, ~, inds] = unique(values);

    % Choose Matlab implementation
    switch mesh.dim
    case 2
        plot_ = @(edges, value, color) plot(...
            reshape(coords(1, edges_to_vertices(:, edges)), 2, []), ...
            reshape(coords(2, edges_to_vertices(:, edges)), 2, []), ...
            'DisplayName', sprintf('%d', value), ...
            'Color', color, 'LineWidth', 3);
    case 3
        plot_ = @(edges, value, color) plot3(...
            reshape(coords(1, edges_to_vertices(:, edges)), 2, []), ...
            reshape(coords(2, edges_to_vertices(:, edges)), 2, []), ...
            reshape(coords(3, edges_to_vertices(:, edges)), 2, []), ...
            'DisplayName', sprintf('%d', value), ...
            'Color', color, 'LineWidth', 3);
    end

    cmap = jet(numel(values));
    plts = gobjects(numel(values), 1);

    % Plot value-by-value
    for j = 1:numel(values)
        edges = marked_edges(inds == j);
        value = values(j);
        color = cmap(j, :);
        p = plot_(edges, value, color);

        % Remember only first object for each value to have concise legend
        plts(j) = p(1);
    end
end


function plts = plot_vertex_markers(mesh, vertex_markers)
    coords = mesh.vertex_coords;

    % Build unique marker values and mapping back to marked vertices
    [marked_vtx, ~, values] = find(vertex_markers(:));
    [values, ~, inds] = unique(values);

    % Choose Matlab implementation
    switch mesh.dim
    case 2
        plot_ = @(verts, value, color) plot(...
            coords(1, verts), coords(2, verts), ...
            '.', 'DisplayName', sprintf('%d', value), ...
            'Color', color, 'MarkerSize', 23);
    case 3
        plot_ = @(verts, value, color) plot3(...
            coords(1, verts), coords(2, verts), coords(3, verts), ...
            '.', 'DisplayName', sprintf('%d', value), ...
            'Color', color, 'MarkerSize', 23);
    end

    cmap = jet(numel(values));
    plts = gobjects(numel(values), 1);

    % Plot value-by-value
    for j = 1:numel(values)
        verts = marked_vtx(inds == j);
        value = values(j);
        color = cmap(j, :);
        p = plot_(verts, value, color);

        % Remember only first object for each value to have concise legend
        plts(j) = p(1);
    end
end
