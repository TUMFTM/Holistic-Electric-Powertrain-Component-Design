function displayBatteryMotorCompatibility(config_motor, battery_GD, secondMotor)
    % Calculate total motor system power
    if secondMotor
        MotorSystemPower = (config_motor.totalPower1 * config_motor.motorsCount1) + ...
                           (config_motor.totalPower2 * config_motor.motorsCount2);
    else
        MotorSystemPower = (config_motor.totalPower1 * config_motor.motorsCount1);
    end    

    % Calculate battery power capability
    batteryPower = battery_GD.SysInfo.U_nom_sys * battery_GD.SysInfo.I_max_sys/1000;

    % Calculate total motor current
    TotalMotorCurrent = config_motor.current1 + config_motor.current2;

    % Get nominal battery voltage
    batteryVoltage = battery_GD.SysInfo.U_nom_sys;

    % Check maximum module voltage compatibility
    if battery_GD.ModInfo.U_max_mod > 60
        disp('Warning: The battery module voltage exceeds 60V.');
    end

    % Improved motor voltage compatibility check: Range-based evaluation
    voltageThresholdLower = 0.9 * batteryVoltage;
    voltageThresholdUpper = 1.1 * batteryVoltage;

    if config_motor.voltage1 < voltageThresholdLower || config_motor.voltage1 > voltageThresholdUpper
        disp(['Warning: Motor 1 voltage (', num2str(config_motor.voltage1), ' V) is outside the optimal range (', ...
              num2str(voltageThresholdLower), ' V - ', num2str(voltageThresholdUpper), ' V) for efficient performance.']);
    end
    if secondMotor
        if config_motor.voltage2 < voltageThresholdLower || config_motor.voltage2 > voltageThresholdUpper
            disp(['Warning: Motor 2 voltage (', num2str(config_motor.voltage2), ' V) is outside the optimal range (', ...
                  num2str(voltageThresholdLower), ' V - ', num2str(voltageThresholdUpper), ' V) for efficient performance.']);
        end
    end
    % Check total motor system power against battery power capability
    if MotorSystemPower > 1.1 * batteryPower
        disp(['Warning: Total motor system power (', num2str(MotorSystemPower), ' kW) exceeds 110% of battery power capability (', ...
              num2str(batteryPower), ' kW).']);
    elseif MotorSystemPower < 0.9 * batteryPower
        disp(['Warning: Total motor system power (', num2str(MotorSystemPower), ' kW) is less than 90% of battery power capability (', ...
              num2str(batteryPower), ' kW).']);
    else
        disp(['Motor system power (', num2str(MotorSystemPower), ' W) is within optimal range relative to battery power capability.']);
    end

    % Compare total power of motors to battery power directly
    totalMotorPower = config_motor.totalPower1 + config_motor.totalPower2;
    if totalMotorPower > batteryPower
        disp(['Warning: Total motor power (', num2str(totalMotorPower), ' W) exceeds the battery power capability (', ...
              num2str(batteryPower), ' W).']);
    end

    % Check total motor current against maximum battery current
    if TotalMotorCurrent > battery_GD.SysInfo.I_max_sys * 1.1
        disp(['Warning: Total motor current (', num2str(TotalMotorCurrent), ' A) exceeds 110% of the battery''s maximum current (', ...
              num2str(battery_GD.SysInfo.I_max_sys), ' A). Risk of overheating.']);
    elseif TotalMotorCurrent < battery_GD.SysInfo.I_max_sys * 0.9
        disp(['Warning: Total motor current (', num2str(TotalMotorCurrent), ' A) is less than 90% of the battery''s maximum current (', ...
              num2str(battery_GD.SysInfo.I_max_sys), ' A). Possible underutilization.']);
    else
        disp(['Total motor current (', num2str(TotalMotorCurrent), ' A) is within optimal range relative to battery current limits.']);
    end
end