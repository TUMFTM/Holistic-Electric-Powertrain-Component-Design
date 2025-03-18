function packComparisonTable = compareBatteryPackData(battery_GD, dataTable)
    % Function to compare battery pack-level attributes between real and replicated data
    %
    % Inputs:
    %   battery_GD - Struct containing replicated battery data
    %   dataTable - Table containing real battery data
    %
    % Output:
    %   packComparisonTable - Table containing the comparison results

    % Initialize comparison results
    vehicleNames = {};
    realPackWeights = [];
    replicatedPackWeights = [];
    realPackVolumes = [];
    replicatedPackVolumes = [];
    realVoltages = [];
    replicatedVoltages = [];
    realCapacities = [];
    replicatedCapacities = [];
    realEnergies = [];
    replicatedEnergies = [];
    weightDeviations = [];
    volumeDeviations = [];
    voltageDeviations = [];
    capacityDeviations = [];
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
            realPackWeights(end+1) = NaN;
            replicatedPackWeights(end+1) = NaN;
            weightDeviations(end+1) = NaN;

            realPackVolumes(end+1) = NaN;
            replicatedPackVolumes(end+1) = NaN;
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

        % Real data: Weight, volume, voltage, capacity, energy
        realPackWeights(end+1) = dataTable{rowIndex, 69}; % Real pack weight
        realPackVolumes(end+1) = dataTable{rowIndex, 80}; % Real pack volume (liters)
        realVoltages(end+1) = dataTable{rowIndex, 70};    % Real pack voltage
        realCapacities(end+1) = dataTable{rowIndex, 72};  % Real pack capacity
        realEnergies(end+1) = dataTable{rowIndex, 74};    % Real pack energy

        % Replicated data: Weight, dimensions, voltage, capacity, energy
        replicatedPackWeights(end+1) = battery_GD(i).battery.SysInfo.mass_sys; % Replicated pack weight
        dim_x_sys = battery_GD(i).battery.SysInfo.dim_x_sys; % Dimension x
        dim_y_sys = battery_GD(i).battery.SysInfo.dim_y_sys; % Dimension y
        dim_z_sys = battery_GD(i).battery.SysInfo.dim_z_sys; % Dimension z

        % Calculate replicated pack volume
        replicatedPackVolumes(end+1) = dim_x_sys * dim_y_sys * dim_z_sys / 1e-3; % Convert to liters
        replicatedVoltages(end+1) = battery_GD(i).battery.SysInfo.U_nom_sys; % Replicated voltage
        replicatedCapacities(end+1) = battery_GD(i).battery.SysInfo.C_sys;   % Replicated capacity
        replicatedEnergies(end+1) = battery_GD(i).battery.SysInfo.E_sys;     % Replicated energy

        % Calculate deviations (only if real value is not NaN or empty)
        if ~isnan(realPackWeights(end))
            weightDeviations(end+1) = 100 * (replicatedPackWeights(end) - realPackWeights(end)) / realPackWeights(end);
        else
            weightDeviations(end+1) = NaN;
        end

        if ~isnan(realPackVolumes(end))
            volumeDeviations(end+1) = 100 * (replicatedPackVolumes(end) - realPackVolumes(end)) / realPackVolumes(end);
        else
            volumeDeviations(end+1) = NaN;
        end

        if ~isnan(realVoltages(end))
            voltageDeviations(end+1) = 100 * (replicatedVoltages(end) - realVoltages(end)) / realVoltages(end);
        else
            voltageDeviations(end+1) = NaN;
        end

        if ~isnan(realCapacities(end))
            capacityDeviations(end+1) = 100 * (replicatedCapacities(end) - realCapacities(end)) / realCapacities(end);
        else
            capacityDeviations(end+1) = NaN;
        end

        if ~isnan(realEnergies(end))
            energyDeviations(end+1) = 100 * (replicatedEnergies(end) - realEnergies(end)) / realEnergies(end);
        else
            energyDeviations(end+1) = NaN;
        end
    end

    % Create a table for the results
    packComparisonTable = table(vehicleNames', realPackWeights', replicatedPackWeights', weightDeviations', ...
        realPackVolumes', replicatedPackVolumes', volumeDeviations', ...
        realVoltages', replicatedVoltages', voltageDeviations', ...
        realCapacities', replicatedCapacities', capacityDeviations', ...
        realEnergies', replicatedEnergies', energyDeviations', ...
        'VariableNames', {'VehicleName', 'RealPackWeight', 'ReplicatedPackWeight', 'WeightDeviation', ...
                          'RealPackVolume', 'ReplicatedPackVolume', 'VolumeDeviation', ...
                          'RealVoltage', 'ReplicatedVoltage', 'VoltageDeviation', ...
                          'RealCapacity', 'ReplicatedCapacity', 'CapacityDeviation', ...
                          'RealEnergy', 'ReplicatedEnergy', 'EnergyDeviation'});

    % Display the comparison table
    disp(packComparisonTable);

    % Optional: Save the comparison table as an Excel file
    writetable(packComparisonTable, 'PackComparisonResults.xlsx');
end