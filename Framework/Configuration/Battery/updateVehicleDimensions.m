function config_vehicle = updateVehicleDimensions(config_vehicle, battery_GD, config_battery)
    % Extract battery scaling factors
    dim_x_factor = config_battery.dim_x_factor;
    dim_y_factor = config_battery.dim_y_factor;
    dim_z_factor = config_battery.dim_z_factor;

    % Convert vehicle dimensions from mm to meters
    vehicle_dim_x_m = config_vehicle.vehicle_dim_x_mm / 1000;
    vehicle_dim_y_m = config_vehicle.vehicle_dim_y_mm / 1000;
    vehicle_dim_z_m = config_vehicle.vehicle_dim_z_mm / 1000;
    
    % Store old frontal area
    old_frontal_area = config_vehicle.frontal_area;
    
    % Initialize frontal area scaling factor
    frontal_area_factor = 1; 

    % Check and update x-dimension
    new_x_m = battery_GD.SysInfo.dim_x_sys / dim_x_factor;
    if vehicle_dim_x_m < new_x_m
        percentage_increase_x = new_x_m / vehicle_dim_x_m;
        config_vehicle.vehicle_dim_x_mm = new_x_m * 1000;
        frontal_area_factor = frontal_area_factor * percentage_increase_x;
        disp('Warning: The vehicle x dimension was increased to fit the battery pack.');
    else
        disp('Vehicle x dimension remains unchanged.');
    end

    % Check and update y-dimension
    new_y_m = battery_GD.SysInfo.dim_y_sys / dim_y_factor;
    if vehicle_dim_y_m < new_y_m
        percentage_increase_y = new_y_m / vehicle_dim_y_m;
        config_vehicle.vehicle_dim_y_mm = new_y_m * 1000;
        frontal_area_factor = frontal_area_factor * percentage_increase_y;
        disp('Warning: The vehicle y dimension was increased to fit the battery pack.');
    else
        disp('Vehicle y dimension remains unchanged.');
    end

    % Apply frontal area scaling only once after x/y updates
    if frontal_area_factor > 1
        config_vehicle.frontal_area = old_frontal_area * frontal_area_factor;
    end

    % Check and update z-dimension
    new_z_m = battery_GD.SysInfo.dim_z_sys / dim_z_factor;
    if vehicle_dim_z_m < new_z_m
        config_vehicle.vehicle_dim_z_mm = new_z_m * 1000;
        disp('Warning: The vehicle z dimension was increased to fit the battery pack.');
    else
        disp('Vehicle z dimension remains unchanged.');
    end
end