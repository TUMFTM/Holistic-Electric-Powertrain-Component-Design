function config_vehicle = updateGearboxWeight(config_vehicle, GB)
    % UPDATEGEARBOXWEIGHT - Updates the weight of the gearbox.
    %
    % Outputs:
    %   - mass_gearbox: Total gearbox weight

    % Access required variables from the base workspace
     config_vehicle.mass_gearbox = GB.m_gearbox;
     config_vehicle.mass_total = config_vehicle.mass_motor + config_vehicle.mass_gearbox + config_vehicle.mass_battery + config_vehicle.mass_rest;

end