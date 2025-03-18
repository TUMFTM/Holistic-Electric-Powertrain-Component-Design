function outputParameters = calculateOutputs(equations, parameters, constants)
    % Initialize output structure
    outputParameters = struct();

    % Iterate through equations to compute output parameters
    for i = 1:length(equations)
        sub_counter = 0; % Counter for substituted variables
        eqn = equations{i}; % Current equation
        
        % Extract left-hand side (LHS) and right-hand side (RHS) of the equation
        lhsEqn = lhs(eqn);
        rhsEqn = rhs(eqn);
        
        % Extract variables present in RHS
        varsInRHS = symvar(rhsEqn);

        % Substitute resolved parameters in the RHS
        for j = 1:length(varsInRHS)
            varName = char(varsInRHS(j)); % Convert symbolic variable to string
            
            if isfield(parameters, varName) && parameters.(varName) ~= 0
                rhsEqn = subs(rhsEqn, sym(varName), parameters.(varName));
                sub_counter = sub_counter + 1;
                
            elseif isfield(outputParameters, varName) && outputParameters.(varName) ~= 0
                rhsEqn = subs(rhsEqn, sym(varName), outputParameters.(varName));
                sub_counter = sub_counter + 1;
                
            elseif isfield(constants, varName) && constants.(varName) ~= 0
                rhsEqn = subs(rhsEqn, sym(varName), constants.(varName));
                sub_counter = sub_counter + 1;
            end
        end
        
        % Solve for output variable only if all required variables are substituted
        if sub_counter == length(varsInRHS)
            outputVar = char(lhsEqn); % Extract LHS variable name
            outputValue = double(rhsEqn); % Compute RHS value
            outputParameters.(outputVar) = outputValue; % Store result
        end
    end    

    % Retrieve constants from workspace
    constants = evalin('base', 'constants');

    %% Special Case: Calculate Maximum Vehicle Speed
    if isfield(outputParameters, 'vehicle_efficiency') && ...
       isfield(parameters, 'motor_peak_power') && ...
       isfield(outputParameters, 'vehicle_mass') && ...
       isfield(constants, 'tire_coefficient') && ...
       isfield(constants, 'density') && ...
       isfield(constants, 'Cw') && ...
       isfield(constants, 'FrontalArea')

        % Assign required constants and parameters
        vehicle_efficiency = outputParameters.vehicle_efficiency;
        motor_peak_power = parameters.motor_peak_power;
        vehicle_mass = outputParameters.vehicle_mass;
        gravitation = constants.gravitation;
        tire_coefficient = constants.tire_coefficient;
        density = constants.density;
        Cw = constants.Cw;
        FrontalArea = constants.FrontalArea;

        % Define resistive power function
        resistivePower = @(v) (vehicle_mass * gravitation * tire_coefficient * v) + ...
                              (0.5 * density * FrontalArea * Cw * v^3);

        % Iteratively estimate max speed
        maxSpeedEstimate = 10; % Initial guess (m/s)
        speedIncrement = 0.1;  % Increment step (m/s)

        while motor_peak_power * vehicle_efficiency > resistivePower(maxSpeedEstimate)
            maxSpeedEstimate = maxSpeedEstimate + speedIncrement;
        end

        % Store maximum speed (converted to km/h)
        outputParameters.vehicle_max_speed = (maxSpeedEstimate - speedIncrement) * 3.6;
    end

    %% Special Case: Calculate Average Energy Consumption
    load('WLTP_class_3.mat');

    if exist('dc', 'var') && isfield(dc, 'speed') && isfield(dc, 'time') && ...
       isfield(outputParameters, 'vehicle_mass') && isfield(outputParameters, 'vehicle_efficiency') && ...
       isfield(constants, 'gravitation') && isfield(constants, 'tire_coefficient') && ...
       isfield(constants, 'density') && isfield(constants, 'FrontalArea') && isfield(constants, 'Cw')

        avg_consumption = calculate_consumption(dc.speed, dc.time, outputParameters.vehicle_mass, ...
                                                constants.gravitation, constants.tire_coefficient, ...
                                                constants.density, constants.FrontalArea, constants.Cw, ...
                                                outputParameters.vehicle_efficiency);

        % Store energy consumption and estimated vehicle range
        outputParameters.vehicle_consumption = avg_consumption;
        outputParameters.vehicle_range = 100 * parameters.battery_energy / (1000 * avg_consumption);
    end

    %% Special Case: Calculate Acceleration Time (0-100 km/h)
    if isfield(outputParameters, 'vehicle_mass') && ...
       isfield(outputParameters, 'vehicle_efficiency') && ...    
       isfield(parameters, 'motor_peak_power') && ...
       isfield(constants, 'gravitation') && ...
       isfield(constants, 'tire_coefficient') && ...
       isfield(constants, 'Cw') && ...
       isfield(constants, 'FrontalArea') && ...
       isfield(constants, 'density')

        % Define parameters
        targetSpeed = 100 / 3.6; % Convert km/h to m/s
        vehicle_mass = outputParameters.vehicle_mass;
        motor_peak_power = parameters.motor_peak_power;
        gravitation = constants.gravitation;
        tire_coefficient = constants.tire_coefficient;
        density = constants.density;
        Cw = constants.Cw;
        FrontalArea = constants.FrontalArea;
        vehicle_efficiency = outputParameters.vehicle_efficiency;

        % Define force functions
        resistiveForce = @(v) vehicle_mass * gravitation * tire_coefficient + ...
                             0.5 * density * FrontalArea * Cw * v^2;

        availableForce = @(v) (motor_peak_power * vehicle_efficiency) / max(v, 0.1); % Avoid division by zero

        % Initialize numerical integration variables
        dt = 0.01; % Time step (seconds)
        speed = 0; % Initial speed
        time = 0;  % Accumulated time

        % Compute acceleration time using numerical integration
        while speed < targetSpeed
            netForce = availableForce(speed) - resistiveForce(speed);
            if netForce <= 0
                break; % Stop if forces balance before reaching target speed
            end
            acceleration = netForce / vehicle_mass;
            speed = speed + acceleration * dt;
            time = time + dt;
        end

        % Store acceleration time result
        if speed >= targetSpeed
            outputParameters.vehicle_acc_time = time;
        else
            outputParameters.vehicle_acc_time = inf; % Acceleration not achievable
        end
    end
end