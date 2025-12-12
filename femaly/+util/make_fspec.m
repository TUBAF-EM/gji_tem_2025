function fs = make_fspec(num, spec)
    % Create string of 'num' comma separated format specifier 'spec'.

    fs = strjoin(repmat({spec}, num, 1), ', ');
end
