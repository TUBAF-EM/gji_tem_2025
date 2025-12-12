function args = intersect_wires_3D(args, tol)
    % Find intersections with wires and add them to wire definitions.
    %
    % Wires are defined as a sequence of straight lines between given
    % points with same ID. Intersection points between the straight line
    % segments of all wires are calculated and added to the point sets of
    % affected wires.
    %
    % This is required to avoid arbitrary renumbering and changing of lines
    % in Gmsh as a result of applying Boolean operations.
    %
    % INPUT PARAMETERS
    %   args    ... Struct, containing geometry information about wires and
    %               domain. Set up in generate_mesh3D.m
    %   tol     ... Scalar, denoting minimum distance for which points are
    %               considered as separate entities.
    %
    % OUTPUT PARAMETERS
    %   args    ... Struct, with updated line definitions.
    %
    % REMARKS
    %   Throughout this function, "wire" will denote the entirety of line
    %   segments of the same ID (defined by two or more points), while
    %   "line" will denote a singular line segment (from one point to one
    %   other).

    % Find intersections with wires and add them to wire definitions.
    if isempty(args.point)
        return;
    end
warning('Intersection needs to be revised!');
% FIXME: revise routine ...  == 2 only works with args.pt!
%        compare with intersetc_wires_2D.m
    is_wire = args.point(:, end) == 2;
    if ~any(is_wire)
        return;
    end

    args = intersect_wire_with_wire(args, is_wire, tol);
    args = intersect_wire_with_point(args, is_wire, tol);
end

function args = intersect_wire_with_wire(args, is_wire, tol)
    % Find wire to wires intersections and add them to wire definitions.
    %
    %   Because the z-coordinate of intersection points at the surface are
    %   calculated within generate_mesh3D.m it is preliminary set to
    %   z = NaN.

    % Get all points of all line entities.
    wi_pt = args.point(is_wire, :);

    % Find number of occurences of line entity point ids.
    [~, id_set] = groupcounts(wi_pt(:, end-1));
    n_wires = length(id_set);

    % Intersect a single wire with all remaining ones.
    for ii = 1:n_wires
        for jj = 1:n_wires-ii
            % Select current wires.
            idx_w1 = wi_pt(:, end-1) == id_set(ii);
            idx_w2 = wi_pt(:, end-1) == id_set(ii+jj);
            wire1 = wi_pt(idx_w1, :);
            wire2 = wi_pt(idx_w2, :);

            if wire1(1, end-2) ~= wire2(1, end-2)
                % Skip intersections between insitu and surface wires case.
                continue;
            end

            if wire1(1, end-2)
                % On surface.
                % Find wire intersections in 2D, i.e. for their projection
                % on x-y plane.
                [p, i ,j] = intersect_w2w(wire1, wire2, args.dim-1, tol);
            else
                % Insitu.
                % Check for intersections in 3D.
                [p, i ,j] = intersect_w2w(wire1, wire2, args.dim, tol);
            end

            if ~isempty(p)
                if wire1(1, end-2)
                    % Set z-coordinate to default.
                    p = [p, NaN*p(:, 1)];
                end

                % Add intersections to point sets of considered wires.
                wire1 = append_point_list(wire1, p, i, id_set(ii));
                wire2 = append_point_list(wire2, p, j, id_set(ii+jj));

                % Append new wire point sets to global point list.
                wi_pt(idx_w1 | idx_w2, :) = [];
                wi_pt = [wi_pt; wire1; wire2];
            else
                continue;
            end
        end
    end

    % Return.
    args.point(is_wire, :) = [];
    args.point = [args.point; wi_pt];
end

function args = intersect_wire_with_point(args, is_wire, tol)
    % Find point to wire intersections and add them to wire definitions.
    %
    % REMARKS
    %   Wires which coincide (in parts) with other wires and therefor are
    %   intersected with wire nodes are also handeled here
    %   (see intersect_wire_with_wire/get_intersection).

    % Get all points of all line entities.
    wi_pt = args.point(is_wire, :);

    % Find number of occurences of line entity point ids.
    [~, id_set] = groupcounts(wi_pt(:, end-1));
    n_wires = length(id_set);

    % Cycle through all straight line segments.
    for ii = 1:n_wires
        % Select current wire.
        idx_w = wi_pt(:, end-1) == id_set(ii);
        wire = wi_pt(idx_w, :);

        % Identify intersection points.
        if wire(1, end-2)   % on surface
            pt_ = args.point(args.point(:, end-2) == 1, :);
            [id_pt, idx_ln] = intersect_w2p(wire(:, 1:args.dim-1), ...
                                               pt_(:, 1:args.dim-1), tol);
        else                % insitu
            pt_ = args.point(args.point(:, end-2) == 0, :);
            [id_pt, idx_ln] = intersect_w2p(wire(:, 1:args.dim), ...
                                               pt_(:, 1:args.dim), tol);
        end

        % Add points to wire point set.
        wire = append_point_list(wire, pt_(id_pt, 1:args.dim), ...
                                         idx_ln, id_set(ii));

        % Append new wire point sets to global point list.
        wi_pt(idx_w, :) = [];
        wi_pt = [wi_pt; wire]; %#ok<*AGROW>
    end

    % Return.
    args.point(is_wire, :) = [];
    args.point = [args.point; wi_pt];
end

%% Helper

function [p, i ,j] = intersect_w2w(wire1, wire2, dim, tol)
    % Calculate intersection points of two wires including inter. at nodes.

    % Loop over wire segments, i.e. straight lines.
    [p, i, j] = deal([]);
    w(1) = warning('off', 'MATLAB:singularMatrix');
    w(2) = warning('off', 'MATLAB:nearlySingularMatrix');
    for ii = 1:size(wire1, 1)-1
        for jj = 1:size(wire2, 1)-1
            % Get line start and end points.
            a = wire1(ii, 1:dim);
            b = wire1(ii+1, 1:dim);
            c = wire2(jj, 1:dim);
            d = wire2(jj+1, 1:dim);

            % Calculate intersection.
            p_ = calculate_intersection(a, b, c, d);
            if ~isempty(p_)
                p = [p; p_]; %#ok<*AGROW>
                i = [i; ii];
                j = [j; jj];
            end
        end
    end
    warning(w);

    % Helper.
    function p = calculate_intersection(a, b, c, d)

        % Sufficient criterion: find intersection for projection on x-y
        % plane.
        m = get_parameter(1);
        if any(isnan(m)) && dim > 2
            % Change to y-z plane in special case of both lines are
            % projected to a single line in x-y plane.
            m = get_parameter(2);
        end

        % Necessary criterion: check if obtained points coincide.
        p1 = a + m(1)*(b - a);
        p2 = c + m(2)*(d - c);
        if norm(p1-p2) < tol && ...    % rounding errors
           ~any(isinf(m)) && ...       % no solution
           ~any(isnan(m))              % lines coincide: case handeled
                                       % in intersect_wire_with_point.m
            p = p1;
        else
            p = [];
        end

        function m = get_parameter(k)
            % Solve a + m(1)*(b - a) =  c + m(2)*(d - c) for a,b,c,d in R^2
            f = (c(k:k+1) - a(k:k+1)).';
            A = [(b(k:k+1) - a(k:k+1)).' -(d(k:k+1) - c(k:k+1)).'];
            m = A\f;
            if any(m > 1) || any(m < 0)
                m = Inf*m; % point doesn't lie between line nodes
            end
        end
    end
end

function [id_pt, i] = intersect_w2p(wire, p, tol)
    % Get distance point to line by exploiting relations within triangle.

    % Loop over straight line segments.
    [id_pt, i] = deal([]);
    for ii = 1:size(wire, 1)-1
        % Get triangle edge lengths.
        a = vecnorm(p             - wire(ii+1, :), 2, 2);
        b = vecnorm(wire(ii+1, :) - wire(ii, :),   2, 2);
        c = vecnorm(p             - wire(ii, :),   2, 2);

        % Get internal angels.
        alpha = acosd((b.^2 + c.^2 - a.^2) ./ 2.*b.*c);
        gamma = acosd((b.^2 + a.^2 - c.^2) ./ 2.*b.*a);

        % Get point to line distance using Heron's formula.
        s = (a+b+c) ./ 2;
        h = 2.*sqrt((s.*(s-a).*(s-b).*(s-c)))./c;

        % Obtain relevant point ids and corresponding line ids.
        id_pt_ = find(real(alpha) < 90 & real(gamma) < 90 & real(h) < tol);
        id_pt = [id_pt; id_pt_];
        i = [i; ii+0*id_pt_];
    end
end

function wire = append_point_list(wire, p, ln_idx, ln_id)
    % Add point to point set with surface/insitu and line id marker.

    dim = size(p, 2);

    % Check if point is already part of wire (intersection at wire node).
    if wire(1, end-2)
        % Note: omit z-coordinate as set to NaN in p.
        is_already_node = ismember(p(:, 1:dim-1), wire(:, 1:dim-1), 'rows');
    else
        is_already_node = ismember(p(:, 1:dim), wire(:, 1:dim), 'rows');
    end
    if all(is_already_node)
        return;
    else
        ln_idx(is_already_node) = [];
        p(is_already_node, :) = [];
    end

    % Add surface/insitu, line id and type markers.
    tmp = 0*p(:, 1);
    if wire(1, end-2)
        p = [p, 1 + tmp, ln_id + tmp, 2 + tmp];
    else
        p = [p,     tmp, ln_id + tmp, 2 + tmp];
    end

    % Sort ascendingly.
    [ln_idx, i] = sort(ln_idx);
    p = p(i, :);

	% Loop over line segment ids.
    i_shift = 0;
    for i = unique(ln_idx).'
        is_ln = ln_idx == i;
        n_p = numel(find(is_ln));
        if n_p > 0
            % Sort points such that they are ordered with increasing
            % distance to first point of line segment.
            p_ = sort_points(p(is_ln, :), i);

        else
            p_ = p(is_ln, :);
        end

        % Add point(s) to wire.
        wire = [wire(1:(i+i_shift), :); p_; wire(i+(1+i_shift):end, :)];

        % Increase shift index by number of just now added points.
        i_shift = i_shift + n_p;
    end

    % Helper.
    function p = sort_points(x, idx)
        % Get distance between points and first node of line segment.
        if wire(1, end-2)
            dist = vecnorm(wire(idx, 1:2) - x(:, 1:2), 2, 2);
        else
            dist = vecnorm(wire(idx, 1:3) - x(:, 1:3), 2, 2);
        end

        % Omit identical points.
        [dist, idx_unique] = unique(dist);
        x = x(idx_unique, :);

        % Sort ascendingly.
        [~, k] = sort(dist);
        p = x(k, :);
    end
end
