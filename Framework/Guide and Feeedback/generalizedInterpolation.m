function predictedValue = generalizedInterpolation(dataTable, independentColumn, dependentColumn, inputValue)
    % Extract relevant columns
    independentData = dataTable.(independentColumn);
    dependentData = dataTable.(dependentColumn);
    
    % Filter rows with non-zero and non-NaN values
    validRows = ~isnan(independentData) & ~isnan(dependentData) & ...
                independentData > 0 & dependentData > 0;
    
    % Extract valid data
    validIndependentData = independentData(validRows);
    validDependentData = dependentData(validRows);
    
    % Handle duplicate values in the independent variable by averaging the corresponding dependent values
    [uniqueIndependentData, ~, idx] = unique(validIndependentData);
    averagedDependentData = accumarray(idx, validDependentData, [], @mean);
    
    % Check if there are enough data points for interpolation
    if length(uniqueIndependentData) < 2
        error('Not enough valid data for interpolation.');
    end
    
    % Create interpolation function
    interpFunc = @(x) interp1(uniqueIndependentData, averagedDependentData, x, 'linear', 'extrap');
    
    % Predict the dependent variable
    predictedValue = interpFunc(inputValue);
    
end