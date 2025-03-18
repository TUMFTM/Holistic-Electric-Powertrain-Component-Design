function avg_consumption = calculate_consumption(dc_speed, dc_time, vehicle_mass, gravitation, tire_coefficient, ...
                                                 density, frontal_area, cw, vehicle_efficiency)
    % Calculates average energy consumption in kWh/100 km
    % Parameters:
    %   dc_speed: array of speeds at each time step (m/s)
    %   dc_time: array of time intervals at each step (s)
    %   vehicle_mass: mass of the vehicle (kg)
    %   gravitation: gravitational constant (m/s^2)
    %   tire_coefficient: coefficient of rolling resistance
    %   density: air density (kg/m^3)
    %   frontal_area: vehicle frontal area (m^2)
    %   cw: drag coefficient
    %   vehicle_efficiency: drivetrain efficiency (0-1)
    %
    % Returns:
    %   avg_consumption: average consumption in kWh/100 km

    % Initialize total energy and distance
    total_energy = 0; % Joules
    total_distance = 0; % Meters

    % Loop through each time step
    for i = 1:length(dc_speed)
        % Get current speed and time
        speed = dc_speed(i); % m/s
        time = dc_time(i); % s

        % Calculate rolling resistance
        rolling_resistance = vehicle_mass * gravitation * tire_coefficient;

        % Calculate aerodynamic drag
        aerodynamic_drag = 0.5 * density * frontal_area * cw * speed^2;

        % Total force
        total_force = rolling_resistance + aerodynamic_drag;

        % Power required (Watts)
        power = (total_force * speed) / vehicle_efficiency;

        % Energy for this step (Joules)
        energy_step = power * time;

        % Add to total energy and distance
        total_energy = total_energy + energy_step;
        total_distance = total_distance + speed * time;
    end

    % Convert total energy to kWh
    total_energy_kwh = total_energy / 3.6e6; % Convert Joules to kWh

    % Convert total distance to kilometers
    total_distance_km = total_distance / 1000; % Meters to kilometers

    % Calculate average consumption in kWh/100 km
    avg_consumption = (total_energy_kwh / total_distance_km) * 100;

end