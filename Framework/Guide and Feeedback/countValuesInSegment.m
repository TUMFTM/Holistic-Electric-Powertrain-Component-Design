function counts = countValuesInSegment(dataTable, segment, columnName, valuesToCount)
    % COUNTVALUESINSEGMENT - Counts occurrences of specific values in a column for a given segment.
    %
    % Inputs:
    %   dataTable     - A MATLAB table containing vehicle data.
    %   segment       - The segment to filter by (e.g., 'A', 'B', 'C').
    %   columnName    - The column name in the table to analyze.
    %   valuesToCount - A cell array of values to count (numeric or string).
    %
    % Outputs:
    %   counts        - A structure containing the counts for each value in valuesToCount.
    
    % Extract the Segment column and the specified column values
    segmentColumn = dataTable{:, 'Segment'}; % Use the column name 'Segment'
    columnValues = dataTable{:, columnName}; % Access the specified column dynamically by name
    
    % Find rows where Segment matches the input
    isSegmentMatch = strcmp(segmentColumn, segment);
    
    % Filter the specified column values for the specified segment
    filteredValues = columnValues(isSegmentMatch);
    
    % Initialize counts as a structure
    counts = struct();
    
    % Check if filteredValues is numeric or cell array of strings
    if isnumeric(filteredValues)
        % Numeric values case
        for i = 1:length(valuesToCount)
            valueToCheck = valuesToCount{i}; % Value to compare
            counts.(genvarname(num2str(valueToCheck))) = sum(filteredValues == valueToCheck); % Count matches
        end
    elseif iscell(filteredValues)
        % String values case
        for i = 1:length(valuesToCount)
            valueToCheck = valuesToCount{i}; % Value to compare
            count = sum(strcmp(filteredValues, valueToCheck)); % Use strcmp for comparison
            counts.(matlab.lang.makeValidName(valueToCheck)) = count; % Store count with a valid field name
        end
    else
        error('Unsupported data type for column "%s". Ensure it is numeric or a cell array of strings.', columnName);
    end
end