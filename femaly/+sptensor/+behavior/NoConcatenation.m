classdef NoConcatenation
    % Basic concatenation definitions for derivation of tensor class.
    %
    % Currently no concatenation are supported.

    methods
        function res = cat(dim, varargin) %#ok<MANU,STOUT,VANUS>
            % TODO: Help.

            error('Concatenation is unsupported.');
        end

        function res = horzcat(varargin) %#ok<STOUT,VANUS>
            % TODO: Help.

            error('Concatenation is unsupported.');
        end

        function res = vertcat(varargin) %#ok<STOUT,VANUS>
            % TODO: Help.

            error('Concatenation is unsupported.');
        end
    end
end
