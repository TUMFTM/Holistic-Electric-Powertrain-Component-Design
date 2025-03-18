function motorWeightComparisonTable = compareMotorWeights(motorTable, dataTable)
    % Initialize result storage
    motorNames = {};
    motorDetails = {}; % New column for storing tokens
    preservedWeights = [];
    estimatedWeights = [];
    deviations = [];
    
    % Define regex patterns for different naming formats
      % Define regex patterns for different motor naming formats
    pattern1 = 'motor_map_.*?_(\d+)Nm_(\d+)kW_\d+Arms'; % Standard format
    pattern2 = 'motor_map_PSM\(\d+\)_(\d+)Nm_(\d+)kW_\d+V_\d+Arms'; % PSM format with parentheses
    pattern3 = 'motor_map_PSM_(\d+)Nm_(\d+)kW_\d+V_\d+Arms_V\d+\(\d+\)'; % New PSM format with voltage and versioning

    % Iterate through all motors in motorTable
    for i = 1:height(motorTable)
        motorName = motorTable{i, "FileName"}; % Extract motor name
        motorNames{end+1} = motorName{1}; % Store the motor name in one column

        % Get preserved weight from function
        preservedWeight = getMechanicalWeight(motorName{1});
        preservedWeights(end+1) = preservedWeight;

        % Extract peak torque and power using regex (try all patterns)
        tokens = regexp(motorName{1}, pattern1, 'tokens');
        if isempty(tokens)
            tokens = regexp(motorName{1}, pattern2, 'tokens');
        end
        if isempty(tokens)
            tokens = regexp(motorName{1}, pattern3, 'tokens');
        end

        if isempty(tokens)
            warning('Invalid motor name format for %s. Skipping...', motorName{1});
            preservedWeights(end) = NaN;
            estimatedWeights(end+1) = NaN;
            deviations(end+1) = NaN;
            motorDetails{end+1} = 'Invalid format'; % Add placeholder
            continue;
        end

        % Extract numerical values
        peakTorque = str2double(tokens{1}{1}); % Peak torque (Nm)
        disp(peakTorque)
        peakPower = str2double(tokens{1}{2});  % Peak power (kW)
        % Determine motor type (PSM or ASM)
        if contains(motorName{1}, 'PSM')
            motorType = 'PSM';
            [frontMotorData, rearMotorData] = extractMotorData(dataTable, 'PSM');
 
        elseif contains(motorName{1}, 'ASM')
            motorType = 'ASM';
            [frontMotorData, rearMotorData] = extractMotorData(dataTable, 'ASM');
        else
            warning('Unknown motor type for %s. Skipping...', motorName{1});
            preservedWeights(end+1) = NaN;
            estimatedWeights(end+1) = NaN;
            deviations(end+1) = NaN;
            motorDetails{end+1} = 'Unknown type'; % Add placeholder
            continue;
        end
        
        % Add motor details to the new column
        motorDetails{end+1} = sprintf('%s, %d Nm, %d kW', motorType, peakTorque, peakPower);

        % Estimate motor weight
        estimatedWeight = estimateMotorWeight(frontMotorData, rearMotorData, peakTorque);
        estimatedWeights(end+1) = estimatedWeight;
        
        % Compute deviation in percentage (only if preservedWeight is valid)
        if ~isnan(preservedWeight) && preservedWeight ~= 0
            deviation = ((estimatedWeight - preservedWeight) / preservedWeight) * 100;
        else
            deviation = NaN;
        end
        deviations(end+1) = deviation;
    end
    
    % Create the output table
    motorWeightComparisonTable = table(motorNames', motorDetails', preservedWeights', ...
        estimatedWeights', deviations', ...
        'VariableNames', {'MotorName', 'Details', 'PreservedWeight', 'EstimatedWeight', 'DeviationPercentage'});
    
    % Display the table
    %disp(motorWeightComparisonTable);
    
    % Optionally, save the table to an Excel file
    writetable(motorWeightComparisonTable, 'MotorWeightComparison.xlsx');

  
   
end


function mechanicalWeight = getMechanicalWeight(motormap_name)
    % GETMECHANICALWEIGHT - Retrieves the mechanical weight for a motor from a data table.
    %
    % Inputs:
    %   - motormap_name: The name of the motor map file (string).
    %
    % Outputs:
    %   - mechanicalWeight: The mechanical weight of the motor (numeric).
    
    % Load motor data table
    try
        loadedData = load('motorDataFull.mat');
        motorTable = loadedData.motorTable;
    catch
        error('Failed to load "motorDataFull.mat". Ensure the file exists and contains "motorTable".');
    end

    % Remove the '.mat' suffix from the input motormap_name
    motormap_name_cleaned = erase(motormap_name, '.mat');

    % Find the row in the table where the FileName matches the cleaned input
    row = strcmp(motorTable.FileName, motormap_name_cleaned);

    % Retrieve the mechanical weight if a match is found
    if any(row)
        mechanicalWeight = motorTable.Thermal(row);
    else
        error('Motor name "%s" not found in the motor table.', motormap_name);
    end
end