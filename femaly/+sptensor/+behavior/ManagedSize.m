classdef ManagedSize
    % Basic size definitions for derivation of tensor class.
    %
    % Currently, most common size methods (available for e.g. matrices in
    % MATLAB) are also available for tensor class.

    methods (Abstract, Access = protected)
        % Queries the size property of the current object.
        %
        % Abstract method - concrete definition is given in Tensor3Base or
        % within one of its subclasses.

        sz = querySize(val)
    end

    methods (Access = protected)
        function n = subscriptedNumel(val, varargin) %#ok<VANUS,MANU,STOUT>
            % TODO: Help.

            error('Subscripted access is unsupported.');
        end
    end

    methods
        function ind = end(val, k, n)
            % TODO: Help.

            sz = querySize(val);
            if k <= length(sz)
                if k < n
                    ind = sz(k);
                else
                    ind = prod(sz(k:end));
                end
            else
                ind = 1;
            end
        end

        function tf = isempty(val)
            % TODO: Help.

            sz = querySize(val);
            tf = any(sz == 0);
        end

        function n = length(val)
            % TODO: Help.

            sz = querySize(val);
            n = max([sz, 0]);
        end

        function n = ndims(val)
            % TODO: Help.

            sz = querySize(val);
            n = max(2, find([0, sz] ~= 1, 1, 'last') - 1);
        end

        function n = numel(val, varargin)
            % TODO: Help.

            if nargin == 1
                sz = querySize(val);
                n = prod(sz);
            else
                n = subscriptedNumel(val, varargin{:});
            end
        end

        function varargout = size(val, dim)
            % TODO: Help.

            sz = querySize(val);
            n = max(2, find([0, sz] ~= 1, 1, 'last') - 1);
            sz = sz(1:n);
            if nargin == 1
                if nargout <= 1
                    varargout = {sz};
                else
                    if nargout < n
                        sz = [sz(1:nargout - 1), prod(sz(nargout:end))];
                    else
                        sz = [sz, ones(1, nargout - n)];
                    end
                    varargout = num2cell(sz);
                end
            else
                if dim <= n
                    varargout = {sz(dim)};
                else
                    varargout = {1};
                end
            end
        end

        function val = subsasgn(val, subs, data) %#ok<INUSD,MANU>
            % TODO: Help.

            error('Subscripted access is unsupported.');
        end

        function varargout = subsref(val, subs) %#ok<INUSD,MANU,STOUT>
            % TODO: Help.

            error('Subscripted access is unsupported.');
        end
    end
end
