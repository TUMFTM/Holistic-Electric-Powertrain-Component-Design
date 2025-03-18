function comparisonTable = compareBatteryData(battery_GD, dataTable)
    % Function to compare cell weight, volume, voltage, capacity, and energy
    % between real and replicated data, and calculate percentage deviations.
    
    % Initialize comparison results
    vehicleNames = {}; 
    realCellWeights = [];
    replicatedCellWeights = [];
    weightDeviations = [];
    realCellVolumes = [];
    replicatedCellVolumes = [];
    volumeDeviations = [];
    realVoltages = [];
    replicatedVoltages = [];
    voltageDeviations = [];
    realCapacities = [];
    replicatedCapacities = [];
    capacityDeviations = [];
    realEnergies = [];
    replicatedEnergies = [];
    energyDeviations = [];

    % Iterate over the vehicles in battery_GD
    for i = 1:length(battery_GD)
        if isempty(battery_GD(i).Vehicle)
            continue; % Skip empty entries
        end

        % Get vehicle name
        vehicleName = battery_GD(i).Vehicle;
        vehicleNames{end+1} = vehicleName;

        % Find the corresponding row in dataTable
        rowIndex = find(strcmp(dataTable{:, 2}, vehicleName));
        if isempty(rowIndex)
            % If no matching row found, add NaN and continue
            realCellWeights(end+1) = NaN;
            replicatedCellWeights(end+1) = NaN;
            weightDeviations(end+1) = NaN;
            realCellVolumes(end+1) = NaN;
            replicatedCellVolumes(end+1) = NaN;
            volumeDeviations(end+1) = NaN;
            realVoltages(end+1) = NaN;
            replicatedVoltages(end+1) = NaN;
            voltageDeviations(end+1) = NaN;
            realCapacities(end+1) = NaN;
            replicatedCapacities(end+1) = NaN;
            capacityDeviations(end+1) = NaN;
            realEnergies(end+1) = NaN;
            replicatedEnergies(end+1) = NaN;
            energyDeviations(end+1) = NaN;
            continue;
        end

        % Real data
        realCellWeights(end+1) = dataTable{rowIndex, 38}; % Real cell weight
        realCellVolumes(end+1) = dataTable{rowIndex, 42}; % Real cell volume (in liters)
        realVoltages(end+1) = dataTable{rowIndex, 35}; % Real cell voltage
        realCapacities(end+1) = dataTable{rowIndex, 36}; % Real cell capacity
        realEnergies(end+1) = dataTable{rowIndex, 37}; % Real cell energy

        % Replicated data
        replicatedCellWeights(end+1) = battery_GD(i).battery.BatPara.physical.m; % Replicated cell weight
        dim_x = battery_GD(i).battery.BatPara.physical.dim_x; % Dimension x
        dim_y = battery_GD(i).battery.BatPara.physical.dim_y; % Dimension y
        dim_z = battery_GD(i).battery.BatPara.physical.dim_z; % Dimension z

        % Calculate replicated cell volume
        if dim_x == dim_y
            % Cylindrical cell volume
            radius = dim_x / 2; % Convert diameter to radius
            replicatedCellVolumes(end+1) = pi * (radius^2) * dim_z * 1e3; % Convert to liters
        else
            % Rectangular prism volume
            replicatedCellVolumes(end+1) = dim_x * dim_y * dim_z * 1e3; % Convert to liters
        end

        replicatedVoltages(end+1) = battery_GD(i).battery.BatPara.electrical.U_nom; % Replicated voltage
        replicatedCapacities(end+1) = battery_GD(i).battery.BatPara.electrical.C_A; % Replicated capacity
        replicatedEnergies(end+1) = battery_GD(i).battery.BatPara.electrical.E; % Replicated energy

        % Calculate percentage deviations
        weightDeviations(end+1) = ((replicatedCellWeights(end) - realCellWeights(end)) / realCellWeights(end)) * 100;
        volumeDeviations(end+1) = ((replicatedCellVolumes(end) - realCellVolumes(end)) / realCellVolumes(end)) * 100;
        voltageDeviations(end+1) = ((replicatedVoltages(end) - realVoltages(end)) / realVoltages(end)) * 100;
        capacityDeviations(end+1) = ((replicatedCapacities(end) - realCapacities(end)) / realCapacities(end)) * 100;
        energyDeviations(end+1) = ((replicatedEnergies(end) - realEnergies(end)) / realEnergies(end)) * 100;
    end

    % Create a table for the results
    comparisonTable = table(vehicleNames', realCellWeights', replicatedCellWeights', weightDeviations', ...
        realCellVolumes', replicatedCellVolumes', volumeDeviations', ...
        realVoltages', replicatedVoltages', voltageDeviations', ...
        realCapacities', replicatedCapacities', capacityDeviations', ...
        realEnergies', replicatedEnergies', energyDeviations', ...
        'VariableNames', {'VehicleName', 'RealCellWeight', 'ReplicatedCellWeight', 'WeightDeviationPercent', ...
                          'RealCellVolume', 'ReplicatedCellVolume', 'VolumeDeviationPercent', ...
                          'RealVoltage', 'ReplicatedVoltage', 'VoltageDeviationPercent', ...
                          'RealCapacity', 'ReplicatedCapacity', 'CapacityDeviationPercent', ...
                          'RealEnergy', 'ReplicatedEnergy', 'EnergyDeviationPercent'});

    % Display the comparison table
    disp(comparisonTable);

    % Optional: Save the comparison table as an Excel file
    writetable(comparisonTable, 'BatteryComparisonResults_Percentage.xlsx');
end