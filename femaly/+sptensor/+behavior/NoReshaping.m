classdef NoReshaping
    % Basic reshaping definitions for derivation of tensor class.
    %
    % Currently no reshapings are supported.

    methods
        function val = ctranspose(val) %#ok<MANU>
            % TODO: Help.

            error('Transpose is unsupported.');
        end

        function val = permute(val, order) %#ok<INUSD,MANU>
            % TODO: Help.

            error('Permutation is unsupported.');
        end

        function val = reshape(val, varargin) %#ok<VANUS,MANU>
            % TODO: Help.

            error('Reshaping is unsupported.');
        end

        function val = transpose(val) %#ok<MANU>
            % TODO: Help.

            error('Transpose is unsupported.');
        end
    end
end
