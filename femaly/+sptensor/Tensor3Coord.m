classdef (Sealed) Tensor3Coord < sptensor.Tensor3Base
    % Sparse rank-3 tensor in coordinate representation.
    %
    % Implements a sparse rank-3 tensor with arbitrary size and nonzero
    % element structure, as it arises in many applications, e.g. when
    % arranging the element mass or stiffness matrices of a finite element
    % discretization in separate slices of the tensor.
    %
    % See also Tensor3Base, Tensor3Coord.

    properties (Access = private)
        size_;
        index_;
        value_;
        cache_;
    end

    methods (Access = protected)
        function s = querySize(val)
            % Implements the abstract method from 'behavior.ManagedSize'.

            s = val.size_;
        end
    end

    % Construction of empty array.
    methods (Static)
        function val = empty(varargin)
            % Constructs an empty sparse rank-3 tensor.
            %
            % SYNTAX
            %
            %   val = sptensor.Tensor3Coord.empty([m[, n[, p[, ...]]]])
            %   val = sptensor.Tensor3Coord.empty(sz)
            %
            % INPUT/OUTPUT PARAMETERS
            %
            %   m   ... Integer scalar giving the number of rows of the
            %           sparse tensor. If this is the only provided input,
            %           also the number of columns of the sparse tensor. If
            %           omitted, defaults to one.
            %   n   ... Integer scalar giving the number of columns of the
            %           sparse tensor. If omitted, defaults to the value
            %           given in 'm'.
            %   p   ... Integer scalar giving the number of slices of the
            %           sparse tensor. If omitted, defaults to one.
            %           Additional inputs following this one can be
            %           specified, but are required to all be one.
            %   sz  ... Non-empty row vector of integers that gives the
            %           number of rows, columns, and (optionally) slices.
            %           More dimensions can be specified, but all those
            %           trailing dimensions need to be one.
            %   val ... Empty sparse tensor of the specified size.

            val = sptensor.Tensor3Coord(sptensor.Tensor3Base.preEmpty(varargin{:}));
        end
    end

    % Construction.
    methods
        function self = Tensor3Coord(sz, varargin)
            % Constructs a sparse rank-3 tensor.
            %
            % SYNTAX
            %
            %   self = Tensor3Coord(sz[, index, value])
            %
            % INPUT/OUTPUT PARAMETERS
            %
            %   sz    ... 3-element row vector of non-negative integers
            %             giving the size of the sparse tensor that is to
            %             be constructed.
            %   index ... n x 3 matrix with every row containing the row,
            %             column, and slice index for the corresponding
            %             nonzero element in 'value'. Note that the same
            %             triplet can appear more than once, which will
            %             result in the summation of the corresponding
            %             valies in 'value', like it is also done by
            %             'accumarray' and 'sparse'. If omitted, yields an
            %             empty sparse tensor.
            %   value ... n x 1 vector with values of the nonzero elements
            %             identified by the triplets in 'index'. If
            %             omitted, yields an empty sparse tensor.
            %   self  ... Resulting sparse tensor.
            %
            % See also accumarray, sparse.

            % Pass-through constructor.
            if nargin == 1 && isa(sz, 'Tensor3Coord')
                self = sz;
                return
            end

            % Process arguments.
            [empty, sz, i, v] = sptensor.Tensor3Base.preConstruct(sz, varargin{:});

            % Process indices and values.
            if ~empty
                % Add values at same index.
                [i, ~, p] = unique(i, 'rows');
                v = accumarray(p, v, [length(p), 1]);

                % Eliminate zeros.
                mask = (v ~= 0);
                i = i(mask, :);
                v = v(mask);
            end

            % Assign properties.
            self.size_ = sz;
            self.index_ = i;
            self.value_ = v;
            self.cache_ = containers.Map();
        end
    end

    % Conversion.
    methods
        function data = struct(self)
            % Converts tensor to struct for compatibility with 'varhash'.
            %
            % SYNTAX
            %
            %   data = struct(self)
            %
            % See also Tensor3Base/struct.

            prop = struct(...
                'size' , self.size_ , ...
                'index', self.index_, ...
                'value', self.value_);
            data = struct(...
                'name', class(self), ...
                'prop', prop       , ...
                'vers', 1          );
        end
    end

    % Tensor-specific functionality.
    methods
        function [index, value] = find(val)
            % Finds indices and values of nonzero elements.
            %
            % SYNTAX
            %
            %   [index, value] = find(val)
            %
            % See also Tensor3Base/find.

            index = val.index_;
            value = val.value_;
        end

        function res = full(val)
            % Converts sparse tensor to full 3D array.
            %
            % SYNTAX
            %
            %   res = full(val)
            %
            % See also Tensor3Base/full.

            % Convert to dense ND array.
            res = accumarray(val.index_, val.value_, val.size_);
        end

        function nz = nnz(val)
            % Returns number of nonzero tensor elements.
            %
            % SYNTAX
            %
            %   nz = nnz(val)
            %
            % See also Tensor3Base/nnz.

            nz = numel(val.value_);
        end

        function s = nonzeros(val)
            % Returns nonzero tensor elements.
            %
            % SYNTAX
            %
            %   s = nonzeros(val)
            %
            % See also Tensor3Base/nonzeros.

            % Return only nonzeros.
            s = val.value_;
        end

        function res = ttm(ten, mat, dim)
            % Computes a tensor-matrix product (ttm = tensor times matrix).
            %
            % SYNTAX
            %
            %   res = ttm(ten, mat, dim)
            %
            % See also Tensor3Base/ttm.

            % Validate tensor and matrix arguments.
            assert(isa(ten, 'sptensor.Tensor3Coord'), ...
                'Expected first argument to be a tensor.');
            assert(isnumeric(mat) && ismatrix(mat), ...
                'Expected second argument to be a matrix.');

            % Check dimension argument and size compatibility.
            assert(isscalar(dim) && any(dim == 1:3), ...
                'Invalid dimension argument.');
            assert(size(mat, 1) == ten.size_(dim), ...
                'Incompatible size of tensor and matrix.');

            % Get size, values, and indices for easier access.
            sz = ten.size_;

            % Get presorted 'i', 'm', and 'v' from cache or generate them.
            dim_key = char('0' + dim);
            if ~isKey(ten.cache_, dim_key)
                dim_order = fliplr([dim, 1:dim - 1, dim + 1:3]);
                [i, p] = sortrows(ten.index_, dim_order);
                m = [true; any(diff(i(:, dim_order(1:2)), 1, 1), 2); true];
                v = ten.value_(p);
                ten.cache_(dim_key) = struct('i', i, 'm', m, 'v', v);
                clear('i', 'm', 'p', 'v');
            end
            dim_ijv = ten.cache_(dim_key);
            ten_i = dim_ijv.i;
            ten_m = dim_ijv.m;
            ten_v = dim_ijv.v;

            % Preprocess matrix.
            mat_sz = size(mat);
            [mat_i, mat_j, mat_v] = find(mat);
            mat_m = [true; logical(diff(mat_j, 1, 1)); true];

            % Result size.
            res_sz = [sz(1:dim - 1), mat_sz(2), sz(dim + 1:3)];
            if size(mat_i, 1) == 0
                res = sptensor.Tensor3Coord(res_sz);
                return
            end

            % Walk index vectors to determine number of nonzeros.
            num_nnz = sptensor.ttmCoordCount(...
                ten_i, ten_m, mat_i, mat_m, dim);

            % Walk index vectors again to populate result.
            [res_i, res_v] = sptensor.ttmCoordApply(...
                ten_i, ten_m, ten_v, mat_i, mat_j, mat_m, mat_v, dim, num_nnz);

            % Construct sparse tensor eliminating duplicate indices.
            res = sptensor.Tensor3Coord(res_sz, res_i, res_v);
        end

        function res = ttv(ten, vec, dim)
            % Computes a tensor-vector product (ttv = tensor times vector).
            %
            % SYNTAX
            %
            %   res = ttv(ten, vec, dim)
            %
            % See also Tensor3Base/ttv.

            % Validate tensor and vector arguments.
            assert(isa(ten, 'sptensor.Tensor3Coord'), ...
                'Expected first argument to be a tensor.');
            assert(isnumeric(vec) && iscolumn(vec), ...
                'Expected second argument to be a column vector.');

            % Check dimension argument and size compatibility.
            assert(isscalar(dim) && any(dim == 1:3), ...
                'Invalid dimension argument.');
            assert(size(vec, 1) == ten.size_(dim), ...
                'Incompatible size of tensor and vector.');

            % Get values and indices for easier access.
            v = ten.value_;
            i = ten.index_;

            % Which two dimensions to keep?
            keep = 1:3;
            keep(dim) = [];

            % Handle sparse and dense vectors differently for performance.
            if issparse(vec)
                % Get 'ten_idx' and 'ten_ptr' from cache or generate them.
                dim_key = sprintf('ttv%d', dim);
                if ~isKey(ten.cache_, dim_key)
                    [ten_ptr, ten_idx] = sort(i(:, dim), 1);
                    ten_ptr = accumarray(ten_ptr, 1, [ten.size_(dim), 1]);
                    ten_ptr = [0; cumsum(ten_ptr, 1)];
                    dim_ttv = struct();
                    dim_ttv.idx = ten_idx;
                    dim_ttv.ptr = ten_ptr;
                    ten.cache_(dim_key) = dim_ttv;
                end
                dim_ttv = ten.cache_(dim_key);
                ten_idx = dim_ttv.idx;
                ten_ptr = dim_ttv.ptr;

                % Dissect given sparse vector into row indices and entries.
                vec_i = find(vec);
                vec_v = nonzeros(vec);

                % Optimize the very common case of a (scaled) unit vector.
                if isscalar(vec_i)
                    % Compute values and restrict indices for nonzero.
                    mask = ten_idx(ten_ptr(vec_i) + 1:ten_ptr(vec_i + 1));
                    v = v(mask) .* vec_v;
                    i = i(mask, keep);
                else
                    % FIXME: Add support for this less common case.
                    % -> ten_idx(ten_ptr(x) + 1:ten_ptr(x + 1))
                    error('Tensor3Coord:ttv:SpSingleNonzeroOnly', [...
                        'Sparse vectors with multiple nonzeros ', ...
                        'are not supported yet.']);
                end
            else
                % Compute new values and restrict indices.
                v = v .* vec(i(:, dim));
                i = i(:, keep);
            end

            % Construct sparse matrix eliminating duplicate indices.
            sz = ten.size_(keep);
            res = sparse(i(:, 1), i(:, 2), v, sz(1), sz(2));
        end
        function res = plus(ten1,ten2)
            assert(isa(ten1, 'sptensor.Tensor3Coord'), ...
                'Expected first argument to be a tensor.');
            assert(isa(ten2, 'sptensor.Tensor3Coord'), ...
                'Expected second argument to be a tensor.');
            assert(all(eq(ten1.size_,ten2.size_)), ...
                'Tensors must be the same size.');

            % Get size, values, and indices for easier access.
            % %%%%%%%%%%%%%%%%%%%%%%%
            sz1 = ten1.size_;
            v1 = ten1.value_;
            i1 = ten1.index_;

            v2 = ten2.value_;
            i2 = ten2.index_;


            index = [i1; i2];
            v = [v1; v2];

            % Construct sparse tensor eliminating duplicate indices and summing them sup.
            res = sptensor.Tensor3Coord(sz1, index, v);
        end

        function res = times(ten, num)
            assert(isa(ten, 'sptensor.Tensor3Coord'), ...
                'Expected first argument to be a tensor.');
            assert(isscalar(num), ...
                'Expected second argument to be a scalar.');

            % Get size, values, and indices for easier access.
            sz = ten.size_;
            v = ten.value_;
            i = ten.index_;

            v = v*num;
            % Construct sparse tensor eliminating duplicate indices.
            res = sptensor.Tensor3Coord(sz, i, v);
        end



    end
end
