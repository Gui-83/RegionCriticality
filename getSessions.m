function validData = getSessions(path)
% Open the file for reading
fileID = fopen(path, 'r');

% Initialize a cell array to store the valid data
validData = {};

while ~feof(fileID)
    line = fgetl(fileID);  % Get the next line
    if isempty(line) || startsWith(line, '%')
        continue;  % Skip empty lines or lines that start with '%'
    end
    % Find the position of the '%' in the line
    idx = find(line == '%', 1);
    if ~isempty(idx)
        line = strtrim(line(1:idx-1));  % Keep only the part before '%'
    end
    
    if endsWith(line, '.xml')
        validData{end+1} = line;  % Store the valid data
    end
end

% Close the file
fclose(fileID);
end

