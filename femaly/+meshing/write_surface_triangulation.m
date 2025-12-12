function [] = write_surface_triangulation(mshf_name, geof_name, args)
    % Convert 2D triangulation from .msh to .geo
    %
    % If topography information is given (args), the z-values of mesh
    % vertices are adjusted using inverse distance weighting approach.
    % If no topography information is given, i.e. args is omitted or all
    % z-values of points in args are equal, heights of vertices in mesh are
    % retained.
    %
    % No physical point markers are set, point entities are only collect
    % in pt_list().
    % Physical line markers are set to it's ID from input (args).

    % Read mesh.
    [mesh, ptm, lnm, mn] = io.read_msh(mshf_name, [0, 1, -1]);
    assert(mesh.dim == 2);
    if nargin < 3 || isempty(args.pt)
        % Set heights == 0.
        vtx = [mesh.vertex_coords.', zeros(mesh.num_entities(0), 1)];
    else
        % Interpolate and adjust vertices z-coordinate using known point
        % heights and topography information.
        vtx = idw(mesh.vertex_coords.', args, 4);
    end

    % Get additional connectivities.
    mesh.init_geometric_queries;
    mesh.compute_connectivity(1, 0);
    ln2vtx = mesh.get_connectivity(1, 0).';
    mesh.compute_connectivity(2, 1);
    cell2ln = mesh.get_connectivity(2, 1).';

    % Get .geo code snippets and summarize.
    geo_code = [write_point(vtx, ptm, mn{1}), newline, ...
                write_line(ln2vtx, lnm), newline, ...
                write_surface(cell2ln)];

    % Write .geo file.
    f_id = fopen(geof_name, 'w');
    fprintf(f_id, geo_code);
    fclose(f_id);
end

%% General helper.

function pt = idw(pt, args, e)
    % Interpolate point z-value by inverse distance weighting.

    % Set exponential weighting factor for distance relations.
    if nargin < 3
        e = 2;
    end

    % Summarize topography information.
    topo = [args.topo; args.pt(args.pt(:, 4) == 1, 1:3)];
    topo(isnan(topo(:, 3)), :) = []; % Remove intersections
    pt = [pt, zeros(size(pt, 1), 1)];
    if isempty(topo) || all(topo(1, 3) == topo(2:end, 3))
        return;
    end

    % Get (horizontal) distances from domain center.
    r = vecnorm(pt(:, 1:2) - args.domain_c(:, 1:2), 2, 2);

    % Get points within domain_ri.
    idx_in_ri = r < args.domain_ri;

    % Get distance relations of points (inside domain_ri) to topo points.
    % Note: If pt coincides with topo, lambda = 1/0 = NaN will be observed.
    lambda = util.pairwise_distance(topo(:, 1:2), pt(idx_in_ri, 1:2)) .^ (-e);

    % Get z-values.
    pt(idx_in_ri, 3) = (lambda.'*topo(:, 3)) ./ sum(lambda)';
    pt_is_topo = isnan(pt(:, 3));

    % Find and replace z=NaN.
    [~, idx_is_topo] = ismember(pt(pt_is_topo, 1:2), topo(:, 1:2), 'rows');
    pt(pt_is_topo, 3) = topo(idx_is_topo, 3);

    % Set height to mean value if outside domain_ri.
    pt(~idx_in_ri, 3) =  mean(topo(:, 3));
end

%% Geocode helper.

function gc = write_point(pt, ptm, mn)
    % Define point entities in .geo syntax.

    % Point definitions.
    gc_pt = cell(size(pt, 1)+1, 1);
    gc_pt{1} = 'pt_id = newp;';
    for i = 1:size(pt, 1)
        gc_pt{i+1} = sprintf(['Point(pt_id+%d) = {', ...
                             util.make_fspec(size(pt, 2), '%.17g'), ...
                             '};'], i-1, pt(i, :));
    end

    % Physical point definition.
    pt_id = unique(ptm(ptm ~= 0));
    is_pt_id = ptm == pt_id.';
    n_pt_id = size(is_pt_id, 2);
    gc_ptm = cell(3*n_pt_id, 1);
    for ii = 1:n_pt_id
        % Prepare physical point list definition (unique coords).
        % -> also define wire/loop corner nodes as physical (dummy) points
        %    to ensure that their ids are preserved after applying
        %    BooleanFragments operations
        cur_id = num2str(pt_id(ii));
        cur_pt_uni = find(is_pt_id(:, ii));
        gc_ptm{3*ii-2} = sprintf(['pt_id_list() = {', ...
                                   util.make_fspec(length(cur_pt_uni), ...
                                                   'pt_id+%d'), ...
                                   '};'], cur_pt_uni-1);
        if strcmp(mn{2}{ii}, 'refine_point')
            % Var 1:
            % Keep line corner/end points defining them as physical point.
            % pt_txt = 'refine_point';
            % gc_ptm{3*ii-1} = ['Physical Point("', pt_txt, ...
                              % '") = pt_id_list();']; % Gmsh assignes id automatically
            % Var 2:
            % Drop traceing line corner/end points, just list them in
            % 'pt_list()' later.
            gc_ptm{3*ii-1} = '';
        else
            pt_txt = ['point_', cur_id];
            gc_ptm{3*ii-1} = ['Physical Point("', pt_txt, ...
                              '", ', cur_id, ') = pt_id_list();'];
        end
        % Append point entity (unique coords) list.
        gc_ptm{3*ii} = 'pt_list() += pt_id_list();';
    end

    % Write geo code.
    gc = strjoin([gc_pt; gc_ptm], newline);
end

function gc = write_line(ln, lnm)
    % Define line entities in .geo syntax.

    % Line definitions.
    gc_ln = cell(size(ln, 1)+1, 1);
    gc_ln{1} = 'ln_id = newl;';
    for i = 1:size(ln, 1)
        id_list = ln(i, :) - 1;
        gc_ln{i+1} = sprintf(['Line(ln_id+%d) = {', ...
                              util.make_fspec(size(ln, 2), 'pt_id+%d'), ...
                              '};'], i-1, id_list);
    end

    % Physical line definition.
    if iscell(lnm)
        gc_lnm = cellfun(@(x) {write_physical_line(x)}, lnm.');
        gc_lnm = vertcat(gc_lnm{:});
    else
        gc_lnm = write_physical_line(lnm);
    end

    % Write geo code.
    gc = strjoin([gc_ln; gc_lnm], newline);
end

function gc = write_surface(tri)
    % Define surface entity in .geo syntax.

    % Line loop definitions.
    n_tri = size(tri, 1);
    gc_lnl = cell(n_tri+1, 1);
    gc_lnl{1} = 'lnl_id = newll;';
    for i = 1:n_tri
        id_list = tri(i, :) - 1;
        gc_lnl{i+1} = sprintf(['Line Loop(lnl_id+%d) = {', ...
                               util.make_fspec(size(tri, 2), 'ln_id+%d'), ...
                               '};'], i-1, id_list);
    end

    % Surface definitions.
    gc_ss = cell(n_tri+1, 1);
    gc_ss{1} = 's_id = news;';
    for i = 1:n_tri
        gc_ss{i+1} = sprintf('Plane Surface(s_id+%d) = {lnl_id+%d};', ...
                             i-1, i-1);
    end

    % Write surface list for physical surface definition.
    gc_sm = sprintf(['bnd_air_earth() = {', ...
                      util.make_fspec(n_tri, 's_id+%d'), '};'], 0:(i-1));

    % Write geo code.
    gc = strjoin([gc_lnl; gc_ss; gc_sm], newline);
end

function gc = write_physical_line(lnm)
    % Write physical line definition.

    if any(lnm)
        % Get markers idx.
        m = full(unique(lnm(lnm ~= 0)));
        idx = lnm == unique(m).';

        % Loop over markers and write physical point definitions.
        gc = cell(numel(m)+1, 1);
        for i = 1:numel(m)
            id_list_ = find(idx(:, i)) - 1;
            gc{i} = sprintf(['Physical Line("line_%d", %d) = {', ...
                             util.make_fspec(numel(id_list_), 'ln_id+%d'), ...
                             '};'], m(i), m(i), id_list_);
        end

        % Write line list.
        id_list_ = find(lnm) - 1;
        gc{end} = sprintf(['ln_list() = {', ...
                  util.make_fspec(numel(id_list_), 'ln_id+%d'), '};'], ...
                  id_list_);
    else
        gc = '';
    end
end
