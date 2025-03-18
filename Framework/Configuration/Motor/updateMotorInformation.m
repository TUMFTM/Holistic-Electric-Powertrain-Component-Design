function [config_motor, config_vehicle] = updateMotorInformation(config_motor, motormap_name_1, motormap_name_2, config_vehicle)
    % UPDATEMOTORINFORMATION - Extracts motor parameters from filenames and updates configurations.
    %
    % Inputs:
    %   config_motor    - Structure storing motor configuration details.
    %   motormap_name_1 - Name of the first motor map file (string).
    %   motormap_name_2 - Name of the second motor map file (string or empty if not applicable).
    %   config_vehicle  - Structure storing vehicle configuration details.
    %
    % Outputs:
    %   config_motor    - Updated motor configuration structure.
    %   config_vehicle  - Updated vehicle configuration structure.

    % Define regex patterns for all three formats
    pattern_old = 'motor_map_.*_(\d+)Nm_(\d+)kW_(\d+)Arms\.mat';
    pattern_new = 'motor_map_(?:\w+_)?(\d+)Nm_(\d+)kW_(\d+)Arms(?:_V\d+\(\d+\))?\.mat';
    pattern_with_parentheses = 'motor_map_(?:\w+\(\d+\)_|\w+_)?(\d+)Nm_(\d+)kW(?:_\d+V)?_(\d+)Arms(?:_V\d+\(\d+\))?\.mat';

    %% Process First Motor
    tokens_1 = regexp(motormap_name_1, pattern_old, 'tokens');
    if isempty(tokens_1)
        tokens_1 = regexp(motormap_name_1, pattern_new, 'tokens'); % Try new pattern if old fails
    end
    if isempty(tokens_1)
        tokens_1 = regexp(motormap_name_1, pattern_with_parentheses, 'tokens'); % Try pattern with parentheses if others fail
    end

    if ~isempty(tokens_1)
        config_motor.torque1   = str2double(tokens_1{1}{1});  % Peak torque (Nm)
        config_motor.power1    = str2double(tokens_1{1}{2});  % Peak power (kW)
        config_motor.current1  = str2double(tokens_1{1}{3});  % Motor current (Arms)

        % Compute total peak values based on the number of motors
        config_motor.totalTorque1 = config_motor.torque1 * config_motor.motorsCount1;
        config_motor.totalPower1  = config_motor.power1 * config_motor.motorsCount1;
    else
        error('Invalid motormap_name_1 format. Expected format: "motor_map_.*_(Torque)Nm_(Power)kW_(Current)Arms.mat".');
    end

    %% Process Second Motor (if applicable)
    if isempty(motormap_name_2)
        % No second motor → Set default values
        config_motor.type2         = '';
        config_motor.torque2       = 0;
        config_motor.power2        = 0;
        config_motor.current2      = 0;
        config_motor.totalTorque2  = 0;
        config_motor.totalPower2   = 0;
    else
        tokens_2 = regexp(motormap_name_2, pattern_old, 'tokens');
        if isempty(tokens_2)
            tokens_2 = regexp(motormap_name_2, pattern_new, 'tokens'); % Try new pattern if old fails
        end
        if isempty(tokens_2)
            tokens_2 = regexp(motormap_name_2, pattern_with_parentheses, 'tokens'); % Try pattern with parentheses if others fail
        end

        if ~isempty(tokens_2)
            config_motor.torque2   = str2double(tokens_2{1}{1}); % Peak torque (Nm)
            config_motor.power2    = str2double(tokens_2{1}{2}); % Peak power (kW)
            config_motor.current2  = str2double(tokens_2{1}{3}); % Motor current (Arms)

            % Compute total peak values
            config_motor.totalTorque2 = config_motor.torque2 * config_motor.motorsCount2;
            config_motor.totalPower2  = config_motor.power2 * config_motor.motorsCount2;
        else
            error('Invalid motormap_name_2 format. Expected format: "motor_map_.*_(Torque)Nm_(Power)kW_(Current)Arms.mat".');
        end
    end
    config_motor.skalieren = false;

    %% Update Vehicle Weight Based on Motor Configurations
    config_vehicle = updateMotorWeight(config_vehicle, config_motor, motormap_name_1, motormap_name_2);
end