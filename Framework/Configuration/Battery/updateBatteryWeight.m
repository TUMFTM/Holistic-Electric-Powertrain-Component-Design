function config_vehicle = updateBatteryWeight(config_vehicle, battery_GD)
    % UPDATEBATTERYWEIGHT - Updates the weight of the battery.
    %
    % Outputs:
    %   - mass_battery: Total battery weight

    % Access required variables from the base workspace

    config_vehicle.mass_battery = battery_GD.SysInfo.mass_sys;
    config_vehicle.mass_total = config_vehicle.mass_motor + config_vehicle.mass_gearbox + config_vehicle.mass_battery + config_vehicle.mass_rest;

end