function solver = get_solver_handle(type, varargin)
    % Creates function handle for fwp solve.
    %
    % INPUT PARAMETER
    %   type ... Char, denoting solver type
    %
    % OPTIONAL PARAMETER
    %   param ... Solver depending additional parameters
    %
    % OUTPUT PARAMETER
    %   solver ... Handle taking matrix A and vector b providing solution
    %              x = A^-1 b

    switch type
        case {'chol', 'backslash'}
            solver = @(A, b) decomposition_solver(A, b, type, varargin);
        case 'mumps'
            solver = @(A, b) mumps_solver(A, b, varargin);
        case 'pcg_amg'
            solver = @(A, b) pcg_amg_solver(A, b, varargin);
        otherwise
            error('Unsupported solver type.');
    end

    function x = mumps_solver(A, b, id)

        % Initialize.
        if ~isempty(id)
            mumps = solving.MUMPS(A, id);
        else
            mumps = solving.MUMPS(A);
        end

        % Solve.
        x = mumps.solve(b);

        % Clean up.
        mumps.delete;
    end

    function x = pcg_amg_solver(A, b, control)

        % Initialize.
        if ~isempty(control)
            amg = solving.HSLMI20(A, control);
        else
            amg = solving.HSLMI20(A);
        end
        x = zeros(size(b));

        % Solve.
        for j = 1:size(b, 2)
            x(:, j) = pcg(A, b(:,j), 1e-10, 100, @amg.solve);
        end

        % Clean up.
        amg.delete;
    end

    function x = decomposition_solver(A, b, type, param)

        % Initialize.
        if strcmp(type, 'chol')
            deco = decomposition(A, type, param{:});
        else
            deco = decomposition(A, 'auto');
        end

        % Solve.
        x = full(deco \ b);

        % Clean up.
        clear('deco');
    end
end
