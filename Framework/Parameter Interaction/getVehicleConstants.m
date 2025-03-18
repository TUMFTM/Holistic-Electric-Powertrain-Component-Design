function constants = getVehicleConstants(vehicle_segment, motor_type)
    
    % Initialize constants structure
    constants = struct();

    % General constants
    
    %constants.const_gearbox = 0.0002;      % Specific volume in m^3/kg         
    %constants.const_acc = 0.25;            % Typical acceleration constant
    constants.gravitation = 9.81;          % Gravitational acceleration (m/s^2)
    constants.tire_coefficient = 0.01;     % Rolling resistance coefficient
    constants.density = 1.2;               % Air density (kg/m³)
    %constants.const_velocity = 30;         % Test velocity (m/s, ~108 km/h)
    constants.stages = 2;                  % Number of gearbox stages
    constants.battery_efficiency = 0.94; 

    % Segment-specific values
    switch vehicle_segment
        case 'A'
            constants.Cw = 0.26;           % Drag coefficient
            constants.FrontalArea = 2.09;   % Frontal area (m²)
            constants.const_mass = 691;            % Base vehicle mass (kg) 
        case 'B'
            constants.Cw = 0.25;
            constants.FrontalArea = 2.19;
            constants.const_mass = 940;            % Base vehicle mass (kg) 
        case 'C'
            constants.Cw = 0.24;
            constants.FrontalArea = 2.36;
            constants.const_mass = 1148;            % Base vehicle mass (kg) 
        case 'D'
            constants.Cw = 0.24;
            constants.FrontalArea = 2.34;
            constants.const_mass = 1286;            % Base vehicle mass (kg) 
        case 'Tesla 3'
            constants.Cw = 0.23;
            constants.FrontalArea = 2.22;
            constants.const_mass = 1306;
        otherwise
            error('Invalid vehicle segment. Choose from Sedan, SUV, Hatchback, Truck.');
    end

    % Motor efficiency based on motor type
    switch motor_type
        case 'PSM'
            %constants.motor_efficiency = 0.95; % Permanent Magnet Synchronous Motor
            constants.motor_specific_weight = 0.25; %      
            constants.motor_corner_speed = 3000;   % typical corner speed in rpm
        case 'ASM'
            %constants.motor_efficiency = 0.90; % Asynchronous Motor
            constants.motor_specific_weight = 0.35; %    
            constants.motor_corner_speed = 2000;   % typical corner speed in rpm
        otherwise
            error('Invalid motor type. Choose from PSM, ASM, SRM.');
    end
    
    % Gearbox efficiency
    %constants.gearbox_efficiency = 0.95; % Fixed as a constant
end