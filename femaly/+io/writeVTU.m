function writeVTU(filename, mesh, varargin)
    % Writes a unstructured grid and some grid quantities to VTK format.
    %
    % SYNTAX
    %
    %   writeVTU(filename, mesh, varargin)
    %
    % INPUT PARAMETER
    %   filemane ... Character denoting path of .vtu file to be stored.
    %                Note: .vtu ending expected
    %   mesh     ... Mesh object from FEMALY
    %
    % OPTIONAL PARAMETER
    %   domain_marker ... Vector [m x 1] of cell-domain-ids
    %   cell_marker   ... Vector [m x 1] of cell parameter values
    %   point_marker  ... Vector [3 x p] of point values
    %                     Note: 3D data is expected
    %
    % OUTPUT
    %   A .vtu file is generated at the desired path which can be
    %   visualized e.g. with PARAVIEW.

    %% Check input.

    % Define possible input keys and its properties checks.
    input_keys = {'domain_marker', 'cell_marker', 'point_marker'};

    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, [], @(x) isvector(x) && ...
                                                    length(x) == size(mesh.cells, 2));
    parser_obj.addParameter(input_keys{2}, [], @(x) isvector(x) && ...
                                                    length(x) == size(mesh.cells, 2));
    parser_obj.addParameter(input_keys{3}, [], @(x) ismatrix(x) && ...
                                                    (size(x, 1) == 1 || size(x, 1) == 3) && ...
                                                    length(x) == size(mesh.vertex_coords, 2));

    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    % Check file extension.
    [~, ~, ext] = fileparts(filename);
    assert(strcmp(ext, '.vtu'), ...
        'Expected file extension ''.vtu'', but ''%s'' found.', ext);

    % Check mesh.
    assert(isa(mesh, 'meshing.Mesh'), 'Expected FEMALY mesh object.');

    %% Fetch mesh info.

    % Set data types.
    [cell_data, point_data] = deal([]);
    if ~isempty(args.cell_marker)
        cell_data.cell_marker = args.cell_marker(:).';
    end
    if ~isempty(args.domain_marker)
        cell_data.domain_marker = args.domain_marker(:).';
    end
    if ~isempty(args.point_marker)
        point_data.point_marker = args.point_marker;
    end

    % Count entities in mesh.
    num_points = size(mesh.vertex_coords, 2);
    num_cells = size(mesh.cells, 2);

    %% Export.

    % Open file.
    vtk = io.VTKXMLWriter(filename, true, 'UnstructuredGrid');
    vtk.openTag('UnstructuredGrid');
    vtk.openTag('Piece', ...
                'NumberOfPoints', uint32(num_points), ...
                'NumberOfCells', uint32(num_cells));

    % Write points.
    vtk.openTag('Points');
    vtk.writeDataArray(mesh.vertex_coords, ...
                      'NumberOfComponents', '3');
    vtk.closeTag('Points');

    % Write cells.
    vtk.openTag('Cells');
    vtk.writeDataArray(int32(mesh.cells - 1), ...
        'Name', 'connectivity');
    vtk.writeDataArray(int32(4 * (1:num_cells)), ...
        'Name', 'offsets');
    vtk.writeDataArray(repmat(uint8(10), 1, num_cells), ...
        'Name', 'types');
    vtk.closeTag('Cells');

    % Write point and cell data.
    vtk.writePointData(point_data, num_points);
    vtk.writeCellData(cell_data, num_cells);

    % Close file.
    vtk.closeTag('Piece');
    vtk.closeTag('UnstructuredGrid');
    vtk.finalize();
end
