classdef MUMPS_OOC < handle

    properties (Access = public)
        id
    end

    properties (GetAccess = public, SetAccess = private)
        mumps
        A
        num
        n_sys
        SYM
        SYM_PERM
    end

    properties (Access = private)
        save_dir
        unique_name
    end

    methods (Access = public)

        function obj = MUMPS_OOC(A, num, id)
            % Factor matrix using MUMPS and store factorization to disc.
            %
            % INPUT ARGS
            %   A   ... sparse square matrix
            %   num ... scalar, file name identifier
            % OPTIONAL ARGS
            %   id ... MUMPS parameters obtained with initmumps()
            % PRIVATE ARGS
            %   SYM      ... symmetry flag, determined from given matrix A
            %   SYM_PERM ... permutation vector, constructed by MUMPS
            %   save_dir ... directory path where mumps factorizations will
            %                be stored in
            %   n_sys    ... size of square sparse matrix
            %   unique_name ... unique object name created from num

            if ~obj.is_installed()
                error('MUMPS package not installed!');
            end

            % Set optimization for symmetric matrices
            % Note: if not exactly symmetric, MUMPS may terminate with
            %       seg. fault!
            if issymmetric(A) && norm(A - A.', 'fro') < 1e-12
                obj.SYM = 1;
            else
                obj.SYM = 0;
            end

            % Define unique file name to ensure that no files will be
            % overwritten
            obj.unique_name = char(java.util.UUID.randomUUID);

            % Set properties
            assert(isscalar(num));
            base = num2str(num);
            obj.num = [base '_' obj.unique_name];
            obj.A = A;
            if isreal(obj.A)
                obj.mumps = @dmumps;
            else
                obj.mumps = @zmumps;
            end

            % Update MUMPS control structure (for given parameters)
            obj.id = initmumps();
            if nargin > 2
                obj.id = update_struct(obj.id, id);
            end

            % Set MUMPS file directory (location for factorization storage)
            save_dir = fullfile(pwd, 'tmp_mumps');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            obj.save_dir = save_dir;
            setenv('MUMPS_SAVE_DIR', obj.save_dir);

            % Set environment variable for object file name
            setenv('MUMPS_SAVE_PREFIX', ['deco_', obj.num]);

            % Initialize MUMPS instance
            obj.id.JOB = -1;
            obj.id = obj.mumps(obj.id);
            obj.id.SYM = obj.SYM;
            obj.check_err();

            % Enable OOC
            obj.id.ICNTL(22) = 1;
            obj.id.ICNTL(35) = 3;

            % Force silent output
            obj.id.ICNTL(1) = -1;
            obj.id.ICNTL(2) = -1;
            obj.id.ICNTL(3) = -1;
            obj.id.ICNTL(4) = 0;

            % Factorize (symbolic and numeric)
            obj.id.JOB = 4; % analyse & factorize
            obj.id = obj.mumps(obj.id, obj.A);
            obj.check_err();

            % Save the instance to disk
            obj.id.JOB = 7; % write to disc
            obj.id = obj.mumps(obj.id, A);

            % Store permutation information
            obj.SYM_PERM = obj.id.SYM_PERM;

            % Replace property with dummy matrix of correct size, required
            % for the obj.mumps()-interface (passing A as second argument
            % is mandatory)
            % Note: original A should not be stored, however, it's
            %       decomposition is written to disc
            obj.n_sys = size(A, 1);
            obj.A = sparse([], [], [], obj.n_sys, obj.n_sys);
        end

        function load(obj)
            % Read factorization from disc.

            % Set environment variable
            % Note: if multiple solver objects are user within a (parfor)
            %       loop, each needs its specific environment variable!
            setenv('MUMPS_SAVE_DIR', obj.save_dir);
            setenv('MUMPS_SAVE_PREFIX', ['deco_', obj.num]);
            if ~isempty(obj.id)
                return
            end

            % Restore flags and permutation
            tmp_id = initmumps();
            tmp_id.JOB = -1;
            tmp_id.N = obj.n_sys;
            tmp_id.ICNTL(22) = 1; % activate MUMPS out-of-core (0 = in RAM)
            tmp_id.ICNTL(35) = 3; % set complete ooc
            obj.id = obj.mumps(tmp_id);

            % Restore saved instance
            obj.id.SYM = obj.SYM;
            obj.id.SYM_PERM = obj.SYM_PERM;
            obj.id.JOB = 8; % load from disc
            obj.id = obj.mumps(obj.id, obj.A);
            obj.check_err();
        end

        function unload(obj)
            % Free memory (without removing files from disc).

            if isempty(obj.id)
                return;
            end
            obj.id.JOB = -2; % free memory
            obj.id = obj.mumps(obj.id, obj.A); % output required
            obj.id = [];
        end

        function delete(obj)
            % Delete stored factorization files from disc.

            try
                % Load instance (set environmet variable)
                obj.load();

                % Remove saved files
                obj.id.JOB = -3; % remove from disc
                obj.id = obj.mumps(obj.id, obj.A);
            catch
                % warning('Files already may have been removed from disc.');
            end

            % Cleanup object
            obj.unload();
        end

        function check_err(obj)
            % Read out MUMPS object err/warn messages and rethrow.

            if isempty(obj.id) || obj.id.INFOG(1) == 0
                return
            end

            msg = sprintf('INFOG(1) = %d\nINFOG = %s', ...
                          obj.id.INFOG(1), mat2str(obj.id.INFOG));

            if obj.id.INFOG(1) > 0
                warning(['MUMPS WARNING: ', msg]);
            elseif obj.id.INFOG(1) < 0
                error(['MUMPS ERROR: ', msg]);
            end
        end

        function [x, id] = solve(obj, b)
            % Load factorization from disc and keep infos in memory.
            %
            % Solve equation A*x = b where b can have multiple columns

            % Set environment variable and if required read factorization
            % from disc
            obj.load();

            % Set RHS and solve
            obj.id.RHS = b;
            obj.id.JOB = 3; % solve
            obj.id = obj.mumps(obj.id, obj.A);
            obj.check_err();

            % Assign return values
            x = obj.id.SOL;
            if nargout > 1
                id = obj.id;
            end
        end

    end

    methods (Static)

        function result = is_installed()
            % Test  MUMPS intialization to check if library is installed.

            try
                [~] = initmumps();
                result = true;
            catch ME
                if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
                    result = false;
                else
                    rethrow(ME);
                end
            end
        end

    end

end


function s = update_struct(s, s_new)
    % Replace solver flags.

%FIXME: May move to +utils

    % Remove fields in s_new from s and raise error if field not in s
    s = rmfield(s, fieldnames(s_new));

    % Concatenate s and s_new
    fields = [fieldnames(s); fieldnames(s_new)];
    s = cell2struct([struct2cell(s); struct2cell(s_new)], fields, 1);
end
