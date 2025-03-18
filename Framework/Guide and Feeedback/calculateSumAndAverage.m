function [sumValue, avgValue, vehicleCount] = calculateSumAndAverage(dataTable, segment, variableName)
    % CALCULATESUMANDAVERAGE - Computes the sum, average, and count of vehicles for a given segment and variable.
    %
    % Inputs:
    %   dataTable    - A MATLAB table containing vehicle data.
    %   segment      - The segment of interest (e.g., 'A', 'B', 'C').
    %   variableName - The name of the column in the table to process.
    %
    % Outputs:
    %   sumValue     - The sum of the variable values for the specified segment.
    %   avgValue     - The average of the variable values for the specified segment.
    %   vehicleCount - The count of non-NaN values in the specified segment.
    
    % Extract the segment and target variable columns
    try
        segmentColumn = dataTable{:, 'Segment'}; % Access the segment column
        targetColumn = dataTable{:, variableName}; % Access the target column dynamically by name
    catch
        error('Invalid column name: Check if "Segment" or the provided variableName exists in the table.');
    end
    
    % Filter rows where the segment matches the input
    isSegmentMatch = strcmp(segmentColumn, segment);
    
    % Extract the corresponding values for the segment
    targetValues = targetColumn(isSegmentMatch);
    
    % Remove NaN values (if any) from the extracted values
    targetValues = targetValues(~isnan(targetValues));
    
    % Calculate sum, average, and count of valid values
    sumValue = sum(targetValues);                % Sum of the values
    avgValue = mean(targetValues);               % Average of the values
    vehicleCount = length(targetValues);         % Number of valid (non-NaN) vehicles
end