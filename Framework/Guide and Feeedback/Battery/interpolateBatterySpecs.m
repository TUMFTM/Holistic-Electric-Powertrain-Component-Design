function interpolateBatterySpecs(dataTable, config_motor, secondMotor)
    % INTERPOLATEBATTERYSPECS - Estimates battery energy and voltage based on motor power.
    %
    % Outputs:
    %   - Console feedback with interpolated battery specifications.

    % Compute total motor power based on configuration
    if secondMotor 
        motorPower = config_motor.power1 * config_motor.motorsCount1 + config_motor.power2 * config_motor.motorsCount2; 
    else
        motorPower = config_motor.power1 * config_motor.motorsCount1;
    end    

    % Interpolate battery specifications based on motor power
    predictedEnergy = generalizedInterpolation(dataTable, 'MaxPower_kW_', 'BatteryPackEnergy_kWh__Overview_', motorPower);
    predictedVoltage = generalizedInterpolation(dataTable, 'MaxPower_kW_', 'BatteryPackNominalVoltage_V__Overview_', motorPower);
    
    % Display interpolated battery specifications
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    disp('🔋  **Interpolated Battery Specifications Based on Motor Power**');
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    fprintf('⚡ **Total Installed Motor Power:** %.2f kW\n', motorPower);
    fprintf('🔋 **Estimated Battery Energy (Interpolated from Database):** %.2f kWh\n', predictedEnergy);
    fprintf('🔌 **Estimated Battery Voltage (Interpolated from Database):** %.2f V\n', predictedVoltage);
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
end