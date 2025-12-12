function [res_i, res_v] = ttmCoordApply(ten_i, ten_m, ten_v, mat_i, mat_j, mat_m, mat_v, dim, num_nnz)
    % Performs a sparse tensor-matrix product.
    %
    % SYNTAX
    %
    %   [res_i, res_v] = ttmCoordApply(ten_i, ten_m, ten_v, mat_i, mat_j, mat_m, mat_v, dim, num_nnz)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   ten_i   ... 3 x n matrix of row, column, and slice indices of all
    %               nonzeros entries of a sparse tensor in coordinate
    %               representation. The rows are presorted such, that the
    %               column specified by 'dim' is alternating fastest.
    %   ten_m   ... 1 x n logical column vector that is synchronized with
    %               'ten_i(:, dim)'. All entries that start a new inner
    %               product are marked with true, while the rest remains
    %               set to false.
    %   ten_v   ... 1 x n numeric column vector that is synchronized with
    %               the rows of 'ten_i' and contains the corresponding
    %               nonzero entries of the sparse tensor.
    %   mat_i   ... 1 x m column vector of row indices of all nonzero
    %               entries of the sparse matrix to multiply with.
    %   mat_j   ... 1 x m column vector of column indices of all nonzero
    %               entries of the sparse matrix to multiply with. Every
    %               entry has a corresponding entry in 'mat_i'.
    %   mat_m   ... 1 x m logical column vector that is synchronized with
    %               'mat_i'. Every entry in 'mat_i' that starts a bew
    %               column of the sparse matrix has true in the
    %               corresponding entry of this vector.
    %   mat_v   ... 1 x m numeric column vector that is synchronized with
    %               'mat_i' and 'mat_j' and contains the corresponding
    %               nonzero entries of the sparse matrix.
    %   dim     ... Scalar integer between 1 and 3 (inclusive) that denotes
    %               the dimension of the tensor along which to perform the
    %               multiplication.
    %   num_nnz ... Scalar integer denoting the number of nonzeros of the
    %               product as returned by 'ttmCoordCount'.
    %   res_i   ... 3 x num_nnz matrix of row, column, and slice indices of
    %               the resulting sparse tensor.
    %   res_v   ... 1 z num_nnz numeric column vector that contains the
    %               nonzeros entries of the result corresponding to the
    %               rows of 'res_i'.
    %
    % REMARKS
    %
    %   This is a helper function that is used internally by 'Tensor3Coord'
    %   to compute the sparse tensor result of a sparse tensor-matrix
    %   product. This is usually preceded by a call to 'ttmCoordCount' to
    %   determine the number of nonzeros of this product.
    %
    %   The idea of this code is to repeatedly walk the index vectors of
    %   the involved tensor and matrix, sorted by the dimension along which
    %   to multiply, i.e. along which to perform the inner products.
    %
    % See also ttmCoordCount, Tensor3Coord.

    res_i = zeros(num_nnz, 3);
    res_v = zeros(num_nnz, 1);
    %
    num_nnz = 0;
    ten_num = size(ten_i, 1);
    ten_cur = 1;
    ten_off = ten_cur;
    mat_num = size(mat_i, 1);
    idx_match = false;
    %
    while ten_cur <= ten_num
        mat_cur = 1;
        mat_off = mat_cur;
        %
        while mat_cur <= mat_num
            if ten_i(ten_cur, dim) == mat_i(mat_cur)
                if idx_match
                    old = res_v(num_nnz, 1);
                else
                    idx_match = true;
                    num_nnz = num_nnz + 1;
                    old = 0;
                    res_i(num_nnz, :) = ten_i(ten_cur, :);
                    res_i(num_nnz, dim) = mat_j(mat_cur);
                end
                res_v(num_nnz, 1) = old + ...
                    ten_v(ten_cur) .* mat_v(mat_cur);
                %
                ten_cur = ten_cur + 1;
                mat_cur = mat_cur + 1;
            elseif ten_i(ten_cur, dim) < mat_i(mat_cur)
                ten_cur = ten_cur + 1;
            else % ten_i(ten_cur, dim) > mat_i(mat_cur)
                mat_cur = mat_cur + 1;
            end
            %
            if mat_m(mat_cur) && mat_cur ~= mat_off
                idx_match = false;
                mat_off = mat_cur;
                ten_cur = ten_off;
            end
            %
            if ten_m(ten_cur) && ten_cur ~= ten_off
                idx_match = false;
                while ~mat_m(mat_cur) || mat_cur == mat_off
                    mat_cur = mat_cur + 1;
                end
                mat_off = mat_cur;
                ten_cur = ten_off;
            end
        end
        %
        idx_match = false;
        while ~ten_m(ten_cur) || ten_cur == ten_off
            ten_cur = ten_cur + 1;
        end
        ten_off = ten_cur;
    end
end
