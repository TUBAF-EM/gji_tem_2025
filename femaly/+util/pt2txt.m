function pt2txt(pt, name)
    % Export point coordinates to .txt

    % Fetch.
    assert(size(pt, 2) == 3);
    assert(ischar(name));

    % Export.
    T = table(pt(:, 1), pt(:, 2), pt(:, 3), 'VariableNames',{'x' 'y' 'z'});
    writetable(T, [name, '.txt'], 'FileType', 'text', 'Delimiter',',');
end

