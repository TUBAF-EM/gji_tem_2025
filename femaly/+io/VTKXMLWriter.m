classdef (Sealed) VTKXMLWriter < handle
    % Helper class for writing data files in VTK's XML format.

    properties (GetAccess = private, SetAccess = private)
        debug_ = false;
        type_map_ = [];
        fid_ = [];
        append_mode_ = true; % <- prefer append mode, inline is slow
        append_data_ = {};
        append_length_ = uint32([]);
        append_offset_ = 0;
        tags_ = {};
    end

    methods
        function self = VTKXMLWriter(filename, append_mode, vtk_type)
            % Opens a XML VTK file for writing.

            % check arguments
            assert(ischar(filename), ...
                'Expected argument ''filename'' to be a character array.');
            assert(isscalar(append_mode) && islogical(append_mode), ...
                'Expected argument ''append_mode'' to be a scalar logical.');
            assert(ischar(vtk_type), ...
                'Expected argument ''vtk_type'' to be a character array.');

            % populate type map
            self.populateTypeMap();

            % open file
            self.fid_ = fopen(filename, 'w+', 'native', 'US-ASCII');
            self.append_mode_ = append_mode;

            % write XML header
            fprintf(self.fid_, '<?xml version="1.0"?>\n');

            % open tag 'VTKFile'
            self.openTag('VTKFile', ...
                         'type', vtk_type, ...
                         'version', '0.1', ...
                         'byte_order', 'LittleEndian');
        end

        function delete(self)
            % Closes the file.

            % close file
            fclose(self.fid_);
            self.fid_ = [];
        end

        function finalize(self)
            % Finalizes the write and closes the file.

            % write append data, if any
            if self.append_mode_ && ~isempty(self.append_data_)
                self.writeAppendedData();
            end

            % close tag
            self.closeTag('VTKFile');

            % close file
            self.delete();
        end

        function openTag(self, tag, varargin)
            % Opens a XML tag and tracks nesting.

            % write an open tag
            self.writeTag(false, tag, varargin);
        end

        function emptyTag(self, tag, varargin)
            % Adds an empty XML tag.

            % write an empty tag
            self.writeTag(true, tag, varargin);
        end

        function closeTag(self, tag)
            % Closes the last open XML tag.

            % pop tag from stack and compare with argument
            assert(strcmp(tag, self.tags_{end}), ...
                'Mismatch between given tag and most recent opened tag.');
            self.tags_(end) = [];

            % write a close tag
            indent = repmat('  ', 1, length(self.tags_));
            fprintf(self.fid_, '%s</%s>\n', indent, tag);
            if self.debug_
                fprintf('%s</%s>\n', indent, tag);
            end
        end

        function writeDataArray(self, data, varargin)
            % Writes a data array element and encodes data.

            % check arguments
            assert(isnumeric(data) && ~issparse(data), ...
                'Expected argument ''data'' to be a dense numeric array.');
            assert(isreal(data), ...
                'Expected argument ''data'' to be a non-complex array.');

            % deduce type string
            if isfield(self.type_map_, class(data))
                type_str = self.type_map_.(class(data));
            else
                error('Unsupported data type.');
            end

            if self.append_mode_
                % compute length
                len = numel(data) * numel(typecast(data(1), 'uint8'));

                % get and update offset (4 = size(uint32) = length header)
                offset = self.append_offset_;
                self.append_offset_ = offset + len + 4;

                % append data
                self.append_data_{end + 1} = data;
                self.append_length_(end + 1) = uint32(len);

                % write tag
                self.emptyTag('DataArray', ...
                              'type', type_str, ...
                              'format', 'appended', ...
                              'offset', uint32(offset), ...
                              varargin{:});
            else
                % convert data and compute length
                data = typecast(data(:), 'uint8');
                len = uint32(numel(data));
                len = typecast(len, 'uint8');

                % encode length and data (base64)
                codec = org.apache.commons.codec.binary.Base64();
                len = reshape(char(codec.encode(len)), 1, []);
                data = reshape(char(codec.encode(data)), 1, []);

                % write tag
                self.openTag('DataArray', ...
                             'type', type_str, ...
                             'format', 'binary', ...
                             varargin{:});
                fprintf(self.fid_, '%s', len);
                fprintf(self.fid_, '%s\n', data);
                self.closeTag('DataArray');
            end
        end

        function writePointData(self, point_data, num_points)
            % Writes point data.

            self.writeDataArrayGroup('PointData', point_data, num_points);
        end

        function writeCellData(self, cell_data, num_cells)
            % Writes cell data.

            self.writeDataArrayGroup('CellData', cell_data, num_cells);
        end
    end

    methods (Access = private)
        function populateTypeMap(self)
            % Populates the MATLAB to VTK type map.

            type_map.single = 'Float32';
            type_map.double = 'Float64';
            type_map.int8 = 'Int8';
            type_map.uint8 = 'UInt8';
            type_map.int16 = 'Int16';
            type_map.uint16 = 'UInt16';
            type_map.int32 = 'Int32';
            type_map.uint32 = 'UInt32';

            self.type_map_ = type_map;
        end

        function writeTag(self, empty_tag, tag, attr)
            % Writes either an opening or empty tag.

            % check arguments
            assert(ischar(tag), ...
                'Expected argument ''tag'' to be a character array.');
            assert(mod(length(attr), 2) == 0, ...
                'Expected additional arguments to be key-value pairs.');
            assert(iscellstr(attr(1:2:end)), ...
                'Expected keys in additional arguments to be character arrays.');

            % construct attributes
            attr_str = '';
            for i = 1:2:length(attr)
                [key, val] = attr{i + (0:1)};
                if ischar(val)
                    attr_str = sprintf('%s %s="%s"', attr_str, key, val);
                elseif isinteger(val) && isscalar(val)
                    attr_str = sprintf('%s %s="%d"', attr_str, key, val);
                elseif isnumeric(val) && isscalar(val)
                    attr_str = sprintf('%s %s="%.15g"', attr_str, key, val);
                else
                    error('Unsupported attribute value type.');
                end
            end

            % suffix
            if empty_tag
                suffix = '/';
            else
                suffix = '';
            end

            % write tag
            indent = repmat('  ', 1, length(self.tags_));
            fprintf(self.fid_, '%s<%s%s%s>\n', indent, tag, attr_str, suffix);
            if self.debug_
                fprintf('%s<%s%s%s>\n', indent, tag, attr_str, suffix);
            end

            % push tag to stack
            if ~empty_tag
                self.tags_{end + 1} = tag;
            end
        end

        function writeAppendedData(self)
            % Flushes appended data at the end of the XML file.

            % write tag and appended data
            self.openTag('AppendedData', ...
                         'encoding', 'raw');
            fprintf(self.fid_, '_');
            for i = 1:length(self.append_data_)
                data_type = class(self.append_data_{i});
                fwrite(self.fid_, self.append_length_(i), 'uint32');
                fwrite(self.fid_, self.append_data_{i}, data_type);
            end
            fprintf(self.fid_, '\n');
            self.closeTag('AppendedData');

            % free memory
            self.append_data_ = {};
            self.append_length_ = uint32([]);
        end

        function writeDataArrayGroup(self, what, data, num)
            % Writes a group of data arrays, e.g. cell/point data.

            % skip section, if no data arrays
            if isempty(data) || isempty(fieldnames(data))
                return
            end

            % process fields
            fields = fieldnames(data);
            self.openTag(what);
            for i = 1:length(fields)
                cur_field = fields{i};
                cur_data = data.(cur_field);
                cur_size = size(cur_data);
                assert(isreal(cur_data) && length(cur_size) == 2, ...
                    'Expected data to be a real numeric matrix.');
                assert(~issparse(cur_data), ...
                    'Expected data to be dense.');
                if all(cur_size == [1, num]) || all(cur_size == [num, 1])
                    % scalar
                    self.writeDataArray(cur_data, ...
                                        'Name', cur_field);
                elseif all(cur_size == [3, num])
                    % vector
                    self.writeDataArray(cur_data, ...
                                        'NumberOfComponents', '3', ...
                                        'Name', cur_field);
                else
                    % unsupported
                    error('Size mismatch. Expected data to be scalar or 3D vector.');
                end
            end
            self.closeTag(what);
        end
    end
end
