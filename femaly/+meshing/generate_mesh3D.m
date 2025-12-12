function [mesh, varargout] = generate_mesh3D(varargin)
    % Stores mesh information of a mesh created by Gmsh.
    %
    % Creates a homogeneous half-space whose surface is formed by
    % incorporating observation points and topography information.
    %
    % Coordinate system: Cartesian right-handed
    %                    i.e. z-axis oriented upwards
    %
    % The center of the domain is defined by
    %   x,y:             args.domain_c(1:2)
    %   if not available (max(t)-min(t))/2, t = x,y-coo of [point; topo]
    %   if not available [0, 0]
    %   z:               args.domain_c(3)
    %   if not available mean(t),           t = z-coo of [point; topo]
    %   if not available 0
    %
    % Supported (tested) Gmsh version: 4.8.3, MSH file format version 2
    %
    % SYNTAX
    %   [mesh, varargout] = generate_mesh3D(varargin)
    %
    % OPTIONAL PARAMETER
    %   point      ... Matrix [m x 3] or [m x 5], with:
    %                  [x, y, z(, is_on_surface, id)]
    %                  Given as [m x 3]: point(s) on surface with id == 1
    %                  Given as [m x 5]: point(s):
    %                                    at air-earth surface:
    %                                       is_on_surface == 1
    %                                    at domain interior:
    %                                       is_on_surface == 0
    %                                    with id(s) taken from 5th column
    %   line       ... Matrix [m x 3] or [m x 5], with:
    %                  [x, y, z(, is_on_surface, id)]
    %                  Given as [m x 3]: line/loop on surface with id == 1
    %                  Given as [m x 5]: line(s):
    %                                    at air-earth surface:
    %                                       is_on_surface == 1
    %                                    at domain interior:
    %                                       is_on_surface == 0
    %                                    with id taken from 5th column
    %                                    describe loop:
    %                                       start/end point occures twice
    %                                    describe wire:
    %                                       no point occures twice
    %   topo       ... Matrix [k x 3] of topography points [x, y, z].
    %   ref        ... Scalar, denoting the number of uniform mesh
    %                  refinement steps.
    %   domain_r   ... Scalar, denoting the radius of domain.
    %   domain_ri  ... Scalar, denoting the (inner) radius < domain_r from
    %                  which all points should have the same z-value
    %                  (mean([point(:, 3); topo(:, 3)])).
    %   dem_r      ... Scalar, denoting radius < domain_ri < domain_r up
    %                  to which elevation information from DEM (topo)
    %                  should be used for inverse distance weighting.
    %                  Note: dem_r > radius of the convex hull of point
    %   domain_c   ... Vector [1 x 3] of center of domain.
    %   size_at_pt ... Scalar, denoting the cell sizes at points.
    %   size_at_wr ... Scalar, denoting the cell sizes at wires.
    %   keep_files ... Logical, denoting if .geo and .msh files shouldn't
    %                  be deleted.
    %   shape      ... Character, denoting domain shape (sphere or cuboid).
    %   keep_air   ... Boolean, denoting if air halfspace should be
    %                  preserved.
    %   marker     ... Vector: -1 <= x_i <= 3, of geometric entities to be
    %                  additionally exported from mesh:
    %                   -1 -> vector of physical entity marker
    %                    0 -> vector of point    entity marker
    %                    1 -> vector of line     entity marker
    %                    2 -> vector of face     entity marker
    %                    3 -> vector of volume   entity marker
    %   geo_code   ... Character of arbitrary additional geo code.
    %
    % OUTPUT PARAMETER
    %   mesh      ... Object from Mesh class.
    %   varargout ... Entity marker vectors, ordered by given input marker
    %                 sorting.
    %
    % REMARKS
    %   First, a squared meshed 2D surface is created by Gmsh that only
    %   incorporates points and lines at the surface (meshing is adapted
    %   to these entities).
    %   Then, heights of all mesh nodes are adjusted using inverse distance
    %   weighting, incorporating additional args.topo information.
    %   The resulting surface in 3D space is intersected with the domain of
    %   shape args.shape and meshed in 3D by Gmsh and insitu point or line
    %   entities are added.
    %   If no topography information is provided, i.e. all z-values of
    %   point are equal and optional parameter topo is omitted, surface is
    %   created at z == 0.

    %% Check input.

    % Define possible input keys and its properties checks.
    input_keys = {'ref', 'point', 'line', 'topo', 'domain_r', 'domain_c', ...
                  'size_at_pt', 'marker', 'keep_files', 'shape', ...
                  'keep_air', 'dem_r', 'domain_ri', 'size_at_wr', ...
                  'geo_code'};
    assertRef = @(x) assert(isscalar(x) && ~islogical(x) && x >= 0);
    assertPos = @(x) assert(isempty(x) || (ismatrix(x) && (size(x, 2) == 3 || ...
        (size(x, 2) == 5 && all((x(:, 4) == 1 | x(:, 4) == 0))) ...
                         && all(x(:, 5) > 0))), ...
        ['Matrix [n x 3] (or [n x 5]) with columns ', ...
        '[x, y, z(, insitu, id)] expected.']);
    assertTopo = @(x) assert(ismatrix(x) && size(x, 2) == 3);
    assertScalar = @(x) isscalar(x) && x > 0;
    assertMarker = @(x) isvector(x) && (all(x >= -1) && all(x <= 3));
    assertShape = @(x) ischar(x) && any(strcmp(x, {'sphere', 'cuboid'}));

    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, 0, assertRef);
    parser_obj.addParameter(input_keys{2}, [], assertPos);
    parser_obj.addParameter(input_keys{3}, [], assertPos);
    parser_obj.addParameter(input_keys{4}, [], assertTopo);
    parser_obj.addParameter(input_keys{5}, 5e3, assertScalar);
    parser_obj.addParameter(input_keys{6}, [], assertPos);
    parser_obj.addParameter(input_keys{7}, 10, assertScalar);
    parser_obj.addParameter(input_keys{8}, [], assertMarker);
    parser_obj.addParameter(input_keys{9}, false, @islogical);
    parser_obj.addParameter(input_keys{10}, 'sphere', assertShape);
    parser_obj.addParameter(input_keys{11}, false, @islogical);
    parser_obj.addParameter(input_keys{12}, [], assertScalar);
    parser_obj.addParameter(input_keys{13}, [], assertScalar);
    parser_obj.addParameter(input_keys{14}, 50, assertScalar);
    parser_obj.addParameter(input_keys{15}, '', @ischar);

    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    % Sanity checks.
    % assert(args.domain_r/args.size_at_pt <= 5e4);
    assert(numel(args.topo) <= 1e6, ['Matrix topo [k x 3] too large, ',...
                                     'reduce k to less than 1e6.']);

    % Prepare paths
    path = fileparts(mfilename('fullpath'));
    tmpfile3D = strcat(path, filesep, 'geocode', filesep, 'mesh3D.in');
    tmpfile_surf2D = strcat(path, filesep, 'geocode', filesep, 'mesh3D_surf2D.in');
    geofile_surf3D = strcat(path, filesep, 'mesh_3D_surf3D.geo');
    geofile2D = strcat(path, filesep, 'mesh_3D_surf2D.geo');
    mshfile2D = strcat(path, filesep, 'mesh_3D_surf2D.msh');
    geofile3D = strcat(path, filesep, 'mesh_3D.geo');
    mshfile3D = strcat(path, filesep, 'mesh_3D.msh');

    %% Fetch.

    % Set supported geometry entities to be included in mesh
    % (see inputParser).
    args.entities = {'point', 'line'};

    % Standardize inter.
    args = standardize_input(args);

    % Set domain center and remove entities lying outside given domain size.
    args = set_domain_center(args);

    % Intersect wires and points.
    intersection_tol = 1e-6; % Chose much higher than expected OpenCascade
                             % accuracy of ~1e-8!
    args = meshing.intersect_wires_3D(args, intersection_tol);

    % Set auxiliary radii.
    args.domain_ri_ratio = 0.5;
    args.dem_r_ratio = 0.9;
    args = set_domain_ri(args);
    args = set_dem_r(args);

    % Write .geo code for domain parameter.
    domain_c_geo = [num2str(args.domain_c(1), '%.17g'), ', ', ...
                    num2str(args.domain_c(2), '%.17g'), ', ', ...
                    num2str(args.domain_c(3), '%.17g')];
    domain_r_geo = num2str(args.domain_r,'%.17g');

    % Build geo code for shape.
    switch args.shape
        case 'sphere'
            shape_geo = "1";
        case 'cuboid'
            shape_geo = "2";
        otherwise
            error('Unsupported shape.');
    end

    %% Create 2D surface mesh.

    if ~isempty(args.pt)
        % Sanity check.
        if all(args.pt(1, 3) == args.pt(2:end, 3)) && ...
           any(args.pt(:, 4)) && ~all(args.pt(:, 4)) && ...
           isempty(args.topo)
            args.pt(:, 4) = 1;
            warning('MESHING:generate_mesh3D:SimilarPointAtSurfAndInsitu', ...
                    ['Given points have same height but some ', ...
                     '"is_on_surface" were falsely set to 0.']);
        end

        % Remove altitude information from entity points.
        pt_surf = args.pt(args.pt(:, 4) == 1, :);
        pt_surf = [pt_surf(:, 1:2), ...
                   zeros(size(pt_surf, 1), 1), ...
                   pt_surf(:, 4:end)];

        % Get unique points (coords) of all entities at surface.
        [~, global2unique, unique2global] = unique(pt_surf(:, 1:3), 'rows');

        % Build geo code for points.
        point_geo = write_point(pt_surf, global2unique, unique2global);

        % Build geo code for wires and loops.
        line_geo = write_line(pt_surf, unique2global);
    else
        [point_geo, line_geo] = deal('');
    end

    % Build geo code for threshold filter and refinement.
    cell_size_at_point_geo = num2str(args.size_at_pt, '%.17g');
    cell_size_at_wire_geo = num2str(args.size_at_wr, '%.17g');
    refine_geo = strjoin(repmat({'RefineMesh;'}, 1, args.ref), newline);

    % Insert geo code parts into template.
    geo_template = fileread(tmpfile_surf2D);
    geo_code2D = sprintf(geo_template, ...
                         point_geo, ...
                         line_geo, ...
                         domain_r_geo, ...
                         domain_c_geo, ...
                         shape_geo, ...
                         cell_size_at_point_geo, ...
                         cell_size_at_wire_geo);

    % Write .geo file.
    f = fopen(geofile2D, 'w');
    fprintf(f, geo_code2D);
    fclose(f);

    % Run gmsh.
    cmd = sprintf('gmsh %s -save -format msh2', geofile2D);
    util.run_sys_cmd(cmd);

    %% Create 3D surface mesh.

    % Write .geo file.
    meshing.write_surface_triangulation(mshfile2D, geofile_surf3D, args);
    if ~args.keep_files
        delete(geofile2D, mshfile2D);
    end

    %% Create 3D domain mesh.

    if ~isempty(args.pt)
        % Get insitu entities.
        pt_insitu = args.pt(args.pt(:, 4) == 0, :);

        % Get unique points of all entities within volume.
        [~, global2unique, unique2global] = unique(pt_insitu(:, 1:3), 'rows');

        % Build geo code for points.
        insitu_point_geo = write_point(pt_insitu, global2unique, unique2global);


        % Build geo code for wires and loops.
        insitu_line_geo = write_line(pt_insitu, unique2global);
    else
        [insitu_point_geo, insitu_line_geo] = deal('');
    end

    % Insert .geo code parts into template.
    geo_template = fileread(tmpfile3D);
    geo_code = sprintf(geo_template, ...
                       num2str(args.domain_r, '%.17g'), ...
                       domain_c_geo, ...
                       shape_geo, ...
                       num2str(args.keep_air, '%.17g'), ...
                       insitu_point_geo, ...
                       insitu_line_geo, ...
                       cell_size_at_point_geo, ...
                       cell_size_at_wire_geo, ...
                       args.geo_code, ...
                       refine_geo);

    % Write .geo file.
    f = fopen(geofile3D, 'w');
    fprintf(f, geo_code);
    fclose(f);

    % Run Gmsh.
    cmd = sprintf('gmsh %s -save -format msh2', geofile3D);
    util.run_sys_cmd(cmd);

    % Read 3D mesh.
    varargout = cell(size(args.marker));
    [mesh, varargout{:}] = io.read_msh(mshfile3D, args.marker);

    % Clean up.
    if ~args.keep_files
        delete(geofile3D, mshfile3D, geofile_surf3D);
    end
end

%% General helper.

function args = standardize_input(args)
    % Ensure format of given geomety entities description to be [m x 6].
    %
    % [x, y, z(, is_on_surface, id), type] with
    % type == 1 - point
    %      == 2 - line
    %
    % Note: id == 0 are not supported in Gmsh.

    % Store uniformized info in separate variable.
    n_ent = length(args.entities);
    pt = cell(n_ent, 1);
    for i = 1:n_ent
        if ~isempty(args.(args.entities{i}))
            n = size(args.(args.entities{i}), 1);
            if size(args.(args.entities{i}), 2) == 3
                % Only point coordinates are given: Set all entity points
                % on surface and id == 1 and add type identifier.
                pt{i} = [args.(args.entities{i}), ones(n, 1), ones(n, 1), ...
                           i+zeros(n, 1)];
            else
                % Just add type identifier.
                assert(all(args.(args.entities{i})(:, 5) > 0));
                pt{i} = [args.(args.entities{i}), i+zeros(n, 1)];
            end
        end

        % Ensure validity of higher dim entities.
        if i > 1 && ~isempty(pt{i})
            id = unique(pt{i}(:, 5));
            is_id = pt{i}(:, 5) == id(:).';
            for j = 1:size(is_id, 2)
                tmp = pt{i}(is_id(:, j), 4);
                % Check entity to consist of min. set of required points.
                assert(length(tmp) >= i, ...
                       ['Dim %d object (id %d) consists of too few ', ...
                        'points'], i, id(j));

                % Check entity to be defined either on surface or insitu.
                assert(all(tmp(1) == tmp(2:end)), ...
                       ['Dim %d object (id %d) with parts on surface ', ...
                        'and domain interior not supported.'], i, id(j));
            end
        end
    end
    args.pt = cell2mat(pt);
    args.dim = 3;
end

function args = set_domain_center(args)
    % Set domain center and remove geometry entities outside domain extent.
    %
    % Points:           only those lying outside are removed
    % n-D (lines, ...): entire entities which are lying outside even in
    %                   pieces are removed

    % If not given, get center from points on surface or set to [0, 0, 0].
    if isempty(args.domain_c)
        if ~isempty(args.pt) && any(args.pt(:, 4) == 1)
            point_at_surf = args.pt(logical(args.pt(:, 4)), 1:3);

            % Get centre of surface describing points.
            args.domain_c = [get_mid(point_at_surf(:, 1)), ...
                             get_mid(point_at_surf(:, 2)), ...
                             mean(point_at_surf(:, 3), 1)];
        else
            % Set to default [0, 0, 0].
            args.domain_c = zeros(1, 3);
        end
    end

    % Define location test.
    switch args.shape
            case 'sphere'
                check = @(pt) vecnorm(pt(:, 1:3) - args.domain_c, 2, 2) > ...
                        args.domain_r;
            case 'cuboid'
                check = @(pt) any(abs(pt(:, 1:3) - args.domain_c) > ...
                            args.domain_r, 2);
    end
    if ~isempty(args.pt)
        % Get all entity points lying outside.
        is_outside = check(args.pt);

        % Remove single points or entire n-D entities.
        n_o = length(args.entities);
        is_nD = args.pt(:, end) == 1:n_o;
        for i = 1:n_o
            if i > 1
                % Get id of higher dim. entities (partially) lying outside.
                id_outside = unique(args.pt(is_outside & is_nD(:, i), 5));
                % Identify all related points of those entities.
                is_outside(is_nD(:, i)) = any(args.pt(is_nD(:, i), 5) == ...
                                              id_outside.', 2);
            end
            if any(is_outside(is_nD(:, i)))
                warning('MESHING:generate_mesh3D:EntityPointOutOffDomain', ...
                        '%s(s) exceeding domain extent omitted.', args.entities{i});

            end
        end
        args.pt(is_outside, :) = [];
    end
    if ~isempty(args.topo)
        % Get all topographic points lying outside.
        is_outside = check(args.topo);
        if any(is_outside)
            warning('MESHING:generate_mesh3D:TopoPointOutOffDomain', ...
                    'Topography point outside domain extent are omited.');
            args.topo(is_outside) = [];
        end
    end

    % Helper.
    function out = get_mid(in)
        % Calculate centerpoint of 'in' from its max/min extent.

        [min_in, max_in] = bounds(in);
        out = mean([min_in, max_in]);
    end
end

function args = set_domain_ri(args)
    % Set dist. from domain center up to which topo smooths out to const. height.

    % Check input or set default.
    if isempty(args.domain_ri)
        args.domain_ri = args.domain_ri_ratio * args.domain_r;
    else
        assert(args.domain_r > args.domain_ri, ['Radius of non-flat ',...
            'domain (domain_ri) needs to be smaller than radius of ', ...
            'full domain (domain_r).']);
    end
end

function args = set_dem_r(args)
    % Set dist. from domain center up to which topo info is considered.

    % Check input or set default.
    if isempty(args.dem_r)
        args.dem_r = args.dem_r_ratio * args.domain_ri;
    else
        assert(args.dem_r < args.domain_ri, ['Radius dem_r should be ',...
            'smaller than radius of non-flat domain (domain_ri).']);
    end

    % Remove all topo points outside of dem_r.
    if ~isempty(args.topo)
        is_inside = vecnorm(args.topo(:, 1:2) - args.domain_c(1:2), 2, 2) ...
                            < args.dem_r;
        if any(~is_inside)
            warning('MESHING:generate_mesh3D:TopoPointOutOffDEMRadius', ...
                    'Topography point(s) outside of dem_r are removed.');
            args.topo(~is_inside, :) = [];
        end
    end
end

%% Geo code helper.

function gc = write_point(inp, g2u, u2g)
    % Defines points and point lists in .geo syntax to .geo file template.
    %
    % The point set is reduced to unique occurences.
    % Physical point markers are set, according to point entity ids.
    % All point entities are additionally collected in pt_list().

    if isempty(inp)
        gc = '';
        return;
    end
    pt_uni = inp(g2u, :); % only unique coord points
    n_pt_unique = size(pt_uni, 1);

    % Create point (unique coords) definitions
    % Note: contains only unique points for all entities of all types.
    tmp_gc = cell(n_pt_unique, 1);
    for i = 1:n_pt_unique
        tmp_gc{i} = sprintf(['Point(pt_id+%d) = {', ...
                               util.make_fspec(3, '%.17g'), ...
                               '};'], i-1, pt_uni(i, 1:3));
    end

    % Get point entities (global coords) from list of all points.
    point_ent_global = find(inp(:, 6) == 1);
    n_point_ent = length(point_ent_global);
    % Check if similar point entities in global coords exist in unique
    % coords.
    tmp = u2g(point_ent_global);
    assert(length(unique(tmp)) == n_point_ent);

    % Append point entity (unique coords) list.
    tmp_gc = [tmp_gc; sprintf(['pt_list() += {', ...
                               util.make_fspec(n_point_ent, 'pt_id+%d'), ...
                               '};'], tmp-1)];

	% Get point entity (global coords) ids.
    id_point_ent = unique(inp(point_ent_global, 5));
    n_point_id = length(id_point_ent);
    is_id_point_ent = inp(point_ent_global, 5) == id_point_ent.';

    % Loop over point entities.
    for ii = 1:n_point_id
        % Prepare physical point list definition (unique coords).
        cur_id = num2str(id_point_ent(ii));
        cur_pt_uni = tmp(is_id_point_ent(:, ii));
        tmp_gc = [tmp_gc; sprintf(['pt_id_list() = {', ...
                                   util.make_fspec(length(cur_pt_uni), ...
                                                   'pt_id+%d'), ...
                                   '};'], cur_pt_uni-1)];
        tmp_gc = [tmp_gc; ['Physical Point("', ['point_', cur_id], ...
                           '", ', cur_id, ') = pt_id_list();']];
    end

    % Append topo/insitu point (includes points from wires/loops) lists.
    is_topo = logical(inp(g2u, 4));
    if any(is_topo)
        % topo
        id_topo = find(is_topo);
        tmp_gc = [tmp_gc; sprintf(['topo_pt_list() = {', ...
                                   util.make_fspec(length(id_topo), ...
                                   'pt_id+%d'), '};'], ...
                                   id_topo-1)];
    end
    if any(~is_topo)
        % insitu
        id_insitu = find(~is_topo);
        n_pt_ents_insitu = length(id_insitu);
        tmp_gc = [tmp_gc; sprintf(['pt_list() += {', ...
                                   util.make_fspec(n_pt_ents_insitu, ...
                                   'pt_id+%d'), '};'], ...
                                   id_insitu-1)];
        % Add point to earth volume.
        tmp_gc = [tmp_gc; ['BooleanFragments{Volume{earth_id}; Delete;}', ...
                           '{Point{pt_list()}; Delete;}']];
    end

    % Write geo code.
    gc = strjoin(tmp_gc, newline);
end

function gc = write_line(inp, u2g)
    % Defines wires/loops and line lists in .geo syntax.
    %
    % Physical line markers are set to it's ID from input (args).
    % The straight line set is reduced to unique occurences.

    gc = '';
    if isempty(inp)
        return;
    end

    % Get line entities (global coords) from list of all points.
    line_ent_global = find(inp(:, 6) == 2);
    if ~any(line_ent_global)
        return;
    end

    % Get line entity (global coords) ids.
    id_line_ent = unique(inp(line_ent_global, 5));
    is_id_line_ent = inp(line_ent_global, 5) == id_line_ent.';
    n_line_ent = size(is_id_line_ent, 2);

    % Get unique point ids for each wire/loop.
    tmp = u2g(line_ent_global);
    tmp = arrayfun(@(x) {tmp(is_id_line_ent(:, x))}, 1:n_line_ent);

    % Get straight line point pairs (unique coords) and set gmsh geometric
    % line entity ids for those.
    ln2pt_uni = cellfun(@(x) {sort([x(1:end-1), x(2:end)], 2)}, tmp);
    ln2pt_uni = vertcat(ln2pt_uni{:});
    ln_gid_global = mat2cell((1:size(ln2pt_uni, 1)).', ...
                            cellfun(@length, tmp).'-1, 1);

    % Get unique geometric line entity ids.
    [ln_unique, ~, ln_global2unique] = unique(ln2pt_uni, 'rows');
    ln_gid_unique = cellfun(@(x) {ln_global2unique(x)}, ln_gid_global);

    % Straight line definitions.
    % Note: Straight line definitions occure only once!
    n_ln_unique = size(ln_unique, 1);
    tmp_gc = cell(n_ln_unique+1, 1);
    tmp_gc{1} = 'ln_id = newl;';
    for i = 1:n_ln_unique
        tmp_gc{i+1} = sprintf('Line(ln_id+%d) = {pt_id+%d, pt_id+%d};', ...
                              i-1, ln_unique(i, 1)-1, ln_unique(i, 2)-1);
    end
    tmp_gc = [tmp_gc; 'ln_list() += {', ...
                       sprintf(util.make_fspec(n_ln_unique, ...
                               'ln_id+%d'), (1:n_ln_unique)-1), ...
                       '};'];

    % Loop over wire/loop entities.
    for ii = 1:n_line_ent
        % Sanity check.
        is_at_surf = logical(inp(is_id_line_ent(:, ii), 4));
        assert(all(is_at_surf) == any(is_at_surf), ...
               ['Line entity parts are defined on surface and at ', ...
                'domain interior.']);

        % Prepare physical line definition.
        cur_id = num2str(id_line_ent(ii));
        tmp_gc = [tmp_gc; 'wire_ln_list() = {', ...
                          sprintf(util.make_fspec(numel(ln_gid_unique{ii}), ...
                          'ln_id+%d'), ln_gid_unique{ii}-1), '};']; %#ok<*AGROW>
        tmp_gc = [tmp_gc; ['Physical Curve("', ['line_', cur_id], ...
                           '", ', cur_id, ') = wire_ln_list();']];

        if ~is_at_surf
            % Add wire/loop to earth volume.
            tmp_gc = [tmp_gc; ['BooleanFragments{Volume{earth_id}; Delete;}', ...
                               '{Curve{wire_ln_list()}; Delete;}']];
        end
    end

    % Write geo code.
    gc = strjoin(tmp_gc, newline);
end
