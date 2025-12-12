function run_taylor_test(solver, m_0, S, verbosity)
    % Check taylor series expansion 1st order convergence.

    if nargin < 4
        verbosity = false;
    end

    % Get starting values.
    d_obs_0 = solver(m_0);

    % Initialize test.
    % rng(0815);
    d_m = randn(size(m_0)) .* m_0;
    h_e_min = log10(1e-4);
    h_e_max = log10(1e-1);
    h = logspace(h_e_max, h_e_min, 2*(h_e_max-h_e_min)+1).';

    % Get solutions for all h.
    d_obs = zeros(length(d_obs_0), length(h));
    if verbosity
        fprintf('Solving %d vwp: ', length(h));
    end
    for hh = 1:length(h)
        if verbosity
            str_verbose = fprintf('%d', hh);
        end
        m = m_0 + h(hh)*d_m;
        d_obs(:, hh) = solver(m);
        if verbosity
            fprintf(repmat('\b', 1, str_verbose));
        end
    end
    fprintf('\n');

    % Calculate norms.
    err_n1 = deal(zeros(length(h), 1));
    for hh = 1:length(h)
        % 1th-order (numeric sensitivity)
        err_n1(hh) = norm(d_obs(:, hh)-(d_obs_0(:) + S*h(hh)*d_m), 2);
    end

    % Plot convergence.
    if ~isnan(err_n1(1))
        loglog(h, err_n1/err_n1(1), 'o-r');
    else
        loglog(h, err_n1, 'o-r');
    end
    hold on
        loglog(h, h/h(1), '.b');
        loglog(h, h.^2/h(1).^2, '.r');
    hold off
    xlabel('h');
    ylim([1e-7, 1e1]);
    xlim([min(h), max(h)]);
    legend('e_1(h)', ...
           'O(h)', 'O(h^2)', ...
           'Location', 'SouthWest');
    set(gca, 'XDir', 'reverse');

    % Check order of accuracy.
    rate = @(x) x(1:end-1)./x(2:end);
    if numel(find(rate(err_n1) >= (h(1)/h(2))^1.9)) < 4 % too weak?
        warning('Taylor-test failed!');
    end
end
