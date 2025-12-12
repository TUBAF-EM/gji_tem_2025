function writePVD(filename, files, timesteps)
    % Writes a .pvd file for a series of .vtu files.
    %
    % INPUT PARAMETER
    %   filename  ... String, name of the .pvd file
    %   files     ... Cell array of .vtu file names
    %   timesteps ... Vector of time step values
    %
    % OUTPUT PARAMETER
    %   <filename>.pvd ... .pvd file summarizing all given .vtu with time
    %                      stamp added
    %
    % Attention: Ensure correct (numerical) sorting of file names.
    %            -> dir(...) provides ASCII dictionary order!

    % Fetch input.
    if nargin < 3
        timesteps = 1:length(files);
    end
    if ~endsWith(filename, '.pvd')
        filename = [filename, '.pvd'];
    end

    % Create .pvd
    fid = fopen(filename, 'w');
    fprintf(fid, '<?xml version="1.0"?>\n');
    fprintf(fid, '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid, '  <Collection>\n');
    for i = 1:length(timesteps)
        fprintf(fid, '    <DataSet timestep="%e" file="%s"/>\n', ...
                timesteps(i), files{i});
    end
    fprintf(fid, '  </Collection>\n');
    fprintf(fid, '</VTKFile>\n');
    fclose(fid);
end
