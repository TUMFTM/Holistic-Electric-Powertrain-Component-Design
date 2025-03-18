function cost = calcMotorCost(config_motor, secondMotor)
    % CALC_MOTOR_COST - Calculates the cost of the motor(s) based on configuration.
    %
    % Inputs:
    %   config_motor - A structure containing motor configuration:
    %       - type1, type2: 'PMSM' or 'ASM' for first and second motor
    %       - power1, power2: Peak power of each motor in kW
    %       - motorsCount1, motorsCount2: Number of motors of each type
    %
    % Outputs:
    %   cost - Total cost of the motor(s)

    % Define price per kW based on motor type
    price_per_kW = struct('PMSM', 10, 'ASM', 8); % Prices for each motor type

    % Calculate cost for motor 1
    if isfield(price_per_kW, config_motor.type1)
        cost1 = config_motor.power1 * price_per_kW.(config_motor.type1) * config_motor.motorsCount1;
    else
        error('Unsupported motor type: %s. Supported types are "PMSM" and "ASM".', config_motor.type1);
    end

    % Calculate cost for motor 2
    if secondMotor
        if isfield(price_per_kW, config_motor.type2)
            cost2 = config_motor.power2 * price_per_kW.(config_motor.type2) * config_motor.motorsCount2;
        else
            error('Unsupported motor type: %s. Supported types are "PMSM" and "ASM".', config_motor.type2);
        end
    
        % Compute total cost
        cost = cost1 + cost2;
    else
        cost = cost1;
    end
end