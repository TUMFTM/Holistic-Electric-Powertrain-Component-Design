function selectedCells = select_cell(cells, databank)
    % Initialize output variable
    selectedCells = struct();
    validRowIdx = 1; % To keep track of valid rows

    % Define available capacities for each format and chemistry
    availableCapacities = struct(...
        'Cylindrical_NMC', [27, 3, 5], ...
        'Pouch_NMC', [48, 64, 80, 89, 104, 109, 112, 130, 145, 182], ...
        'Prismatic_LFP', [119, 134, 166, 187, 249, 267, 374, 401, 561], ...
        'Prismatic_NMC', [102, 145, 163, 203, 218, 229, 327, 457, 686] ...
    );
     availableCapacities = struct(...
    'Cylindrical_NMC111', [27, 3, 5], ...
    'Cylindrical_NMC622', [27, 3, 5], ...
    'Cylindrical_NMC721', [27, 3, 5], ...
    'Cylindrical_NMC811', [27, 3, 5], ...
    'Pouch_NMC111', [104, 109, 112, 130, 145, 182, 48, 64, 67, 78, 80, 89], ...
    'Pouch_NMC622', [104, 109, 112, 130, 145, 182, 48, 64, 67, 78, 80, 89], ...
    'Pouch_NMC721', [104, 109, 112, 130, 145, 182, 48, 64, 67, 78, 80, 89], ...
    'Pouch_NMC811', [104, 109, 112, 130, 145, 182, 48, 64, 67, 78, 80, 89], ...
    'Prismatic_LFP', [119, 134, 166, 178, 187, 249, 267, 374, 401, 561, 59, 83], ...
    'Prismatic_NMC111', [102, 145, 163, 203, 218, 229, 305, 327, 457, 490, 686, 73], ...
    'Prismatic_NMC622', [107, 153, 172, 214, 229, 241, 321, 344, 481, 516, 722, 76], ...
    'Prismatic_NMC721', [112, 160, 180, 225, 241, 253, 337, 361, 505, 541, 758, 80], ...
    'Prismatic_NMC811', [118, 168, 189, 235, 252, 265, 353, 378, 529, 567, 794, 84] ...
);
    % Iterate through each row in the databank
    for i = 1:height(databank)
        % Check if column 30 (Cell Format) is empty
        if isempty(databank{i, 30}{1})
            continue; % Skip this row
        end

        % Extract required information
        vehicleName = databank{i, 2}{1}; % Vehicle name
        cellFormat = databank{i, 30}{1}; % Cell format
        cellChemistry = databank{i, 31}{1}; % Chemistry (NMC or LFP)
        secondaryChemistry = databank{i, 32}{1}; % Secondary Chemistry
        materialRatio = databank{i, 33}{1}; % Material Ratio (e.g., NMC622)
        cellCapacity = databank{i, 36}; % Cell capacity (actual value)
         %disp(cellFormat)
       % Standardize cell format names
        
        % Skip if chemistry is not NMC or LFP
        if ~ismember(cellChemistry, {'NMC', 'LFP'})
            continue;
        end

        % Handle secondary chemistry for logging
        if ~isempty(secondaryChemistry) && ~strcmp(secondaryChemistry, '')
            selectedCells(validRowIdx).secondaryChemistry = secondaryChemistry;
        end
        

           

        if strcmp(cellChemistry, 'NMC')
            % Extract the material ratio index (e.g., 'NMC622' -> '622')
            %ratioIndex = extract_ratio_index(materialRatio);
            % Append the ratio to the key
             ratioIndex = extract_ratio_index(materialRatio);
             ratioIndex = regexp(materialRatio, '\d{3}(?:\d{3})?', 'match', 'once');
             if isempty(ratioIndex)
                 continue
             end    
            if isnan(ratioIndex) %|| length(num2str(ratioIndex)) == 2
                continue; % Skip this iteration
            end
             %disp(vehicleName)
             %disp(cellChemistry)
            availableRatios = [111, 622, 721, 811];
            [closestRatio, isExactRatio] = find_closest_ratio(ratioIndex, availableRatios);
            % Check if the ratioIndex matches exactly one of the available ratios
 
            % Save variable to indicate if closest ratio was not exact
            selectedCells(validRowIdx).isExactRatio = isExactRatio;
            selectedCells(validRowIdx).closestRatio = closestRatio;
            key = strcat(cellFormat, '_', cellChemistry, num2str(closestRatio));
            
        else
            % For LFP, the key remains the same as before
            key = strcat(cellFormat, '_', cellChemistry);
            selectedCells(validRowIdx).isExactRatio = true;
        end
        
        % Check if the key exists in the availableCapacities structure
        if ~isfield(availableCapacities, key)
            continue; % Skip if no available capacities for this format/chemistry/ratio
        end
        %disp(key)
        % Access the available capacities for the selected key
        availableCapacitiesForCell = availableCapacities.(key);
         %disp(availableCapacitiesForCell) 
        % Select the closest capacity
        [closestCapacity, ~] = find_closest_value(cellCapacity, availableCapacitiesForCell);
       cellFormat = strtrim(string(cellFormat)); % Ensure it's a string and trim spaces
        if strcmpi(cellFormat, 'Cylindrical')
            cellFormat = 'Cyl';
        elseif strcmpi(cellFormat, 'Prismatic')
            cellFormat = 'Pris';
        elseif strcmpi(cellFormat, 'Pouch')
            cellFormat = 'Pouch';
        else
            continue; % Skip unrecognized formats
        end
        % For NMC, also match the ratio
        if strcmp(cellChemistry, 'NMC')
            % Extract ratio index from materialRatio (e.g., 'NMC622' -> 622)
         %   ratioIndex = extract_ratio_index(materialRatio);
         %   availableRatios = [111, 622, 721, 811];
         %   [closestRatio, isExactRatio] = find_closest_value(ratioIndex, availableRatios);

            % Save variable to indicate if closest ratio was not exact
          %  selectedCells(validRowIdx).isExactRatio = isExactRatio;
          %  selectedCells(validRowIdx).closestRatio = closestRatio;

            % Recreate the cell name for NMC
            selectedCellName = strcat('dummy_', cellFormat, '_NMC', num2str(closestRatio), '_', num2str(closestCapacity), 'Ah');
        else
            % Recreate the cell name for LFP
            selectedCellName = strcat('dummy_', cellFormat, '_LFP_', num2str(closestCapacity), 'Ah');
        end


        CapacitypercentageDeviation = 100 * (closestCapacity - cellCapacity) / cellCapacity;
        % Store the selected cell data
        selectedCells(validRowIdx).vehicleName = vehicleName; % Include vehicle name
        selectedCells(validRowIdx).cellFormat = cellFormat;
        selectedCells(validRowIdx).cellChemistry = cellChemistry;
        selectedCells(validRowIdx).selectedCapacity = closestCapacity;
        selectedCells(validRowIdx).selectedCellName = selectedCellName;
        selectedCells(validRowIdx).capacitydeviation = CapacitypercentageDeviation;
        
        % Increment valid row index
        validRowIdx = validRowIdx + 1;
    end
end

function ratioIndex = extract_ratio_index(materialRatio)
    % Extract the numeric index from materialRatio (e.g., 'NMC622' -> 622)
    ratioIndex = str2double(regexp(materialRatio, '\d+', 'match', 'once'));
end

function [closestValue, isExact] = find_closest_value(value, availableValues)
    % Find the closest value and whether it is exact
    [~, idx] = min(abs(availableValues - value));
    closestValue = availableValues(idx);
    isExact = (closestValue == value);
end

function [closestRatio, isExact] = find_closest_ratio(ratio, availableRatios)
    % Convert the target ratio into percentage values
    targetPercentage = convert_to_percentage(ratio);

    % Initialize variables
    minDeviation = Inf;
    closestRatio = [];
    isExact = false;

    % Convert all available ratios into percentage values
    availablePercentages = cellfun(@convert_to_percentage, num2cell(availableRatios), 'UniformOutput', false);

    % Iterate through the available ratios
    for i = 1:length(availablePercentages)
        % Calculate the deviation for this available ratio
        %disp(size(targetPercentage))
        %disp(size(availablePercentages{i}))
        %disp(targetPercentage)
        %disp(availablePercentages{i})
        deviation = sum(abs(targetPercentage - availablePercentages{i}));

        % Update the closest match if deviation is smaller
        if deviation < minDeviation
            minDeviation = deviation;
            closestRatio = availableRatios(i);
            isExact = (deviation == 0); % Exact match if no deviation
        end
    end
end

function percentageVector = convert_to_percentage(ratio)
    % Convert a ratio into a percentage vector
    ratioString = num2str(ratio); % Convert to string
    if length(ratioString) == 3
        % 3-digit format: Each digit is a proportion (e.g., '811' -> [8, 1, 1])
        components = str2double(regexp(ratioString, '\d', 'match'));
        percentageVector = 100 * components / sum(components); % Normalize to percentages
    elseif length(ratioString) == 6
        % 6-digit format: Each pair represents a percentage (e.g., '665672' -> [66, 56, 72])
        components = str2double({ratioString(1:2), ratioString(3:4), ratioString(5:6)});
        percentageVector = components; % Already percentages
    else
        error('Invalid ratio format: %s. Must be 3-digit or 6-digit.', ratioString);
    end
end