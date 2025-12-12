classdef Tensor3Base ...
        < sptensor.behavior.NoConcatenation ...
        & sptensor.behavior.NoReshaping ...
        & sptensor.behavior.ManagedSize
    % Common interface and functionality of sparse rank-3 tensors.
    %
    % Sparse rank-3 tensors are not natively supported by MATLAB, but often
    % arise in the context of FD/FE discretizations where one wishes to
    % retain easy access to the element matrices that arise during system
    % matrix assembly, either because the constant per-element parameter
    % changes often and one wishes to quickly reassemble the system matrix
    % with a new set of parameters or because one wishes to implement the
    % derivation with respect to those parameters, in which case access to
    % the parameter-dependent parts of the system matrix is needed.
    %
    % Note that this class is an abstract base class and cannot be used
    % directly. At least two specializations exist:
    %
    %   * Tensor3Coord is a rather generic implementation which stores the
    %     sparse tensor in coordinate format and supports rank-3 tensors of
    %     arbitrary size.
    %
    % The two most common operations that can be performed on a sparse
    % tensor are the tensor-matrix product ('ttm', literally tensor times
    % matrix) and the tensor-vector product ('ttv', literally tensor times
    % vector). Their signature is as follows:
    %
    %   res = ttm(ten, mat, dim)
    %   res = ttv(ten, vec, dim)
    %
    % The former yields a sparse tensor while the latter yields a sparse
    % matrix. Both calls require the additional parameter 'dim', which can
    % be 1, 2, or 3 and denotes the dimension along which to multiply the
    % given tensor.
    %
    % See also Tensor3Coord, Tensor3Base/ttm, Tensor3Base/ttv.
    %
    % COPYRIGHT
    %   The class was developed and implemented by Martin Afanasjew.

    % Tensor-specific functionality (abstract).
    methods (Abstract)
        % Finds indices and values of nonzero elements.
        %
        % SYNTAX
        %
        %   [index, value] = find(val)
        %
        % INPUT/OUTPUT PARAMETERS
        %
        %   val   ... Sparse tensor that should be decomposed into its
        %             nonzero entries and corresponding indices.
        %   index ... n x 3 matrix with every row containing the row,
        %             column, and slice index matching the corresponding
        %             nonzero element in 'value'.
        %   value ... n x 1 vector with all nonzero elements of the tensor.
        %
        % See also Tensor3Base/nonzeros.
        [index, value] = find(val)

        % Converts sparse tensor to full 3D array.
        %
        % SYNTAX
        %
        %   res = full(val)
        %
        % INPUT/OUTPUT PARAMETERS
        %
        %   val ... Sparse tensor.
        %   res ... Dense 3D array with the same size and contents like the
        %           given sparse tensor.
        %
        % REMARKS
        %
        %   Note that this function is mostly useful for debugging when the
        %   given tensor is very small. Typical sparse tensors are usually
        %   too large to be represented by a dense 3D array.
        %
        % See also full.
        res = full(val)

        % Returns number of nonzero tensor elements.
        %
        % SYNTAX
        %
        %   nz = nnz(val)
        %
        % INPUT/OUTPUT PARAMETERS
        %
        %   val ... Sparse tensor.
        %   nz  ... Integer scalar giving the number of nonzero elements.
        %
        % See also nnz.
        nz = nnz(val)

        % Returns nonzero tensor elements.
        %
        % SYNTAX
        %
        %   s = nonzeros(val)
        %
        % INPUT/OUTPUT PARAMETERS
        %
        %   val ... Sparse tensor.
        %   s   ... Column vector with all nonzero elements of the tensor.
        %
        % See also Tensor3Base/find.
        s = nonzeros(val)

        % Converts tensor to struct for compatibility with 'varhash'.
        %
        % SYNTAX
        %
        %   data = struct(self)
        %
        % INPUT/OUTPUT PARAMETERS
        %
        %   self ... Sparse tensor.
        %   data ... Struct with version information, class name, and
        %            information that uniquely describes the sparse tensor
        %            using only native MATLAB data types.
        %
        % See also varhash.
        data = struct(self)

        % Computes a tensor-matrix product (ttm = tensor times matrix).
        %
        % SYNTAX
        %
        %   res = ttm(ten, mat, dim)
        %
        % INPUT/OUTPUT PARAMETERS
        %
        %   ten ... Sparse tensor.
        %   mat ... Sparse matrix to multiply with.
        %   dim ... Scalar integer denoting the dimension along which to
        %           multiply. One of 1, 2, or 3. Some implementations do
        %           not support all of the dimensions.
        %   res ... Result of the multiplication, which is again a sparse
        %           tensor.
        %
        % See also Tensor3Base/ttv.
        res = ttm(ten, mat, dim)

        % Computes a tensor-vector product (ttv = tensor times vector).
        %
        % SYNTAX
        %
        %   res = ttv(ten, vec, dim)
        %
        % INPUT/OUTPUT PARAMETERS
        %
        %   ten ... Sparse tensor.
        %   vec ... Column vector to multiply with.
        %   dim ... Scalar integer denoting the dimension along which to
        %           multiply. One of 1, 2, or 3.
        %   res ... Result of the multiplication, which is a sparse matrix.
        %
        % See also Tensor3Base/ttm.
        res = ttv(ten, vec, dim)
    end

    % Helpers for construction.
    methods (Static, Access = protected)
        function sz = preEmpty(varargin)
            % Preprocesses arguments for construction of empty tensor.
            %
            % SYNTAX
            %
            %   sz = preEmpty([varargin])
            %
            % INPUT/OUTPUT PARAMETERS
            %
            %   varargin ... Any arguments that are also accepted by
            %                'zeros', 'ones', and related functions.
            %   sz       ... 3-element row vector derived from 'varargin'
            %                that can be passed as the first argument of
            %                class deriving from this one. An error is
            %                reported if the given arguments cannot be
            %                converted to a valid size of a rank-3 tensor.
            %
            % REMARKS
            %
            %   This is intended as a helper for use in derived classes,
            %   and as such is not exposed to the general public. Look at
            %   existing implementations for details on how to use this.

            % Construct size vector depending on argument count.
            switch nargin
                case 0
                    sz = [1, 1];
                case 1
                    sz = varargin{1};
                    assert(isrow(sz), ...
                        'Single size argument must be a row vector.');
                    if isempty(sz)
                        sz = [0, 0];
                    elseif isscalar(sz)
                        sz = [sz, sz];
                    else
                        % Keep size vector 'sz' as is.
                    end
                otherwise
                    assert(all(cellfun(@isscalar, varargin)), ...
                        'Multiple size arguments must be scalars.');
                    sz = [varargin{:}];
            end

            % Pad size vector with singletons.
            sz = [sz, ones(1, max(0, 3 - length(sz)))];

            % Validate size specification.
            assert(all(fix(sz) == sz) && all(sz >= 0), ...
                'Expected size to be given as non-negative integers.');
            assert(all(sz(4:end) == 1), ...
                'Expected all trailing dimensions starting with the fourth to be one.');
            assert(any(sz(1:3) == 0), ...
                'Expected at least one dimension to be zero.');

            % Return size vector as expected by constructor.
            sz = sz(1:3);
        end

        function [empty, sz, i, v] = preConstruct(sz, i, v)
            % Preprocesses arguments for tensor construction.
            %
            % SYNTAX
            %
            %   [empty, sz, i, v] = preConstruct(sz[, i, v])
            %
            % INPUT/OUTPUT PARAMETERS
            %
            %   sz    ... 3-element row vector with the numer of rows,
            %             columns, and slices of the sparse rank-3 tensor.
            %             This argument is not altered, but is checked and
            %             used to check the other arguments. Alternative,
            %             this can be an object of a class derived from
            %             'Tensor3Base' on input, in which case it will be
            %             decomposed into its parts with 'size' and 'find'
            %             and returned via 'sz', 'i', and 'v'.
            %   i     ... If given, an n x 3 matrix with every row
            %             specifying the row, column, and slice index of
            %             the corresponding entry in 'v'.
            %   v     ... If given, an n x 1 vector with the nonzero
            %             entries of the sparse tensor. The corresponding
            %             rows in 'i' specify the position of the entries
            %             in the tensor.
            %   empty ... Logical scalar that is set to true if the given
            %             arguments suggest an empty tensor, or false
            %             otherwise.
            %
            % REMARKS
            %
            %   This is intended as a helper for use in derived classes,
            %   and as such is not exposed to the general public. Look at
            %   existing implementations for details on how to use this.

            % Convert from a compatible class, if sole argument.
            if nargin == 1 && isa(sz, 'sptensor.Tensor3Base')
                [sz1, sz2, sz3] = size(sz);
                [i, v] = find(sz);
                empty = isempty(i) && isempty(v);
                sz = [sz1, sz2, sz3];
                return
            end

            % Validate size and general stuff.
            assert(nargin == 1 || nargin == 3, ...
                'Expected one or three arguments.');
            assert(isrow(sz) && length(sz) == 3, ...
                'Expected size argument to be a row vector of length three.');
            assert(all(fix(sz) == sz) && all(sz >= 0), ...
                'Expected size argument to be a vector of non-negative integers.');

            % Helpful quantity.
            empty = (nargin == 1) || (isempty(i) && isempty(v));

            % Validate indices and values.
            if empty
                i = zeros(0, 3);
                v = zeros(0, 1);
            else
                assert(ismatrix(i) && size(i, 2) == 3, ...
                    'Expected ''index'' to be a matrix with three columns.');
                assert(iscolumn(v), ...
                    'Expected ''value'' to be a column vector.');
                assert(size(i, 1) == size(v, 1), ...
                    'Expected ''index'' and ''value'' to have the same number of rows.');
                assert(all(min(i, [], 1) >= 1), ...
                    'Some indices are out of range (zero or negative).');
                assert(all(max(i, [], 1) <= sz), ...
                    'Some indices are out of range (larger than size).');
            end
        end
    end

    % Output.
    methods
        function disp(val, name)
            % Displays the sparse tensor, without printing its name.
            %
            % SYNTAX
            %
            %   disp(val[, name])
            %
            % INPUT/OUTPUT PARAMETERS
            %
            %   val  ... Sparse tensor that should be printed. Since
            %            printing of a typical sparse tensor is not very
            %            helpful, only some meta data is output. If you
            %            really want to see the entries of the tensor, you
            %            can use 'full' to convert it to MATLAB's dense
            %            tensor format and then print it.
            %   name ... String with the variable name that should be
            %            prefixed to the output. This is for internal use.
            %            You should use 'display' instead, if you want to
            %            include the variable name in your output.
            %
            % See also disp, Tensor3Base/display.

            % Prepare variable name printing.
            if nargin >= 2
                % Handle 'display' case.
                if isempty(name)
                    name = 'ans';
                else
                    % Keep 'name' as is.
                end
            else
                % Handle 'disp' case.
                name = '';
            end

            % Print variable name.
            if ~isempty(name)
                fprintf('%s =\n', name);
            end

            % Print size and class.
            [m, n, p] = size(val);
            fprintf('   Sparse tensor: %d-by-%d-by-%d', m, n, p);
            fprintf(' (nnz = %s,', num2str(nnz(val)));
            fprintf(' class ''%s'')\n', class(val));
        end

        function display(val)
            % Prints the tensor in expressions without trailing semicolon.
            %
            % SYNTAX
            %
            %   display(val)
            %
            % INPUT/OUTPUT PARAMETERS
            %
            %   val ... Sparse tensor that should be printed. Since
            %           printing of a typical sparse tensor is not very
            %           helpful, only some meta data is output. If you
            %           really want to see the entries of the tensor, you
            %           can use 'full' to convert it to MATLAB's dense
            %           tensor format and then print it.
            %
            % See also display, Tensor3Base/disp.

            disp(val, inputname(1));
        end
    end
end
