function J = assemble_perturbation_sensitivity(solver, m_0, pert)

    % Solve initial fwp.
    nm = length(m_0);
%     assert(nm < 2000, 'You wont want to drink that much coffee ...');
    d0 = solver(m_0);

    % Allocate J.
    J = zeros(length(d0), nm);

    % Loop over cells.
    fprintf('Assemble J\n');
    fprintf(sprintf('Loop over nm = %d parameter: ', nm));
    wb = waitbar(0, 'ii = 0', 'Name', 'Run perturbation method');
    for ii = 1:nm
        % Reset parameter.
        m = m_0;

        % Disturb a single model parameter.
        m(ii) = m_0(ii) * (1 + pert);

        % Get disturbed data.
        d = solver(m);

        % sensitivity calculation
        J(:,ii) = (d - d0) ./ (m_0(ii) * pert);

        % Update waitbar and message
        waitbar(ii/nm, wb, sprintf('ii = %d', ii));
    end
    fprintf(' done.\n');
    close(wb);
end


