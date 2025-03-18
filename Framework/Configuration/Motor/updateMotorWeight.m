function config_vehicle = updateMotorWeight(config_vehicle, config_motor, motormap_name_1, motormap_name_2)
    % UPDATEMOTORWEIGHT - Updates the weight of the motors based on configuration.
    %
    % Inputs:
    %   - config_motor: Struct containing motor configuration details
    %   - motormap_name_1: Name of the first motor map file (string).
    %   - motormap_name_2: Name of the second motor map file (string or empty if not applicable).
    %
    % Outputs:
    %   - config_vehicle: Updated vehicle configuration structure with total motor weight

    % Load mechanical weight from the first motor map
    newWeight1 = getMechanicalWeight(motormap_name_1);

    % Calculate the first motor's total weight
    mass_motor1 = newWeight1  * config_motor.motorsCount1 ;

    % Initialize values for the second motor
    mass_motor2 = 0;

    if ~isempty(motormap_name_2)
        % Load mechanical weight from the second motor map
        newWeight2 = getMechanicalWeight(motormap_name_2);

        mass_motor2 = newWeight2 * config_motor.motorsCount2 ;
    end

    % Calculate total motor weight
    total_mass_motor = mass_motor1 + mass_motor2;

    % Update vehicle configuration
    config_vehicle.mass_motor = total_mass_motor;
    config_vehicle.mass_total = config_vehicle.mass_motor + config_vehicle.mass_gearbox + config_vehicle.mass_battery + config_vehicle.mass_rest;
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