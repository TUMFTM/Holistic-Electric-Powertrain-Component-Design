function import_dc_cycle_v2_1(cycle_name)
    config_vehicle = evalin('base', 'config_vehicle');
    switch cycle_name
        case 'custom'
            % dc = create_acc_cycle(config_vehicle.zero_to_corner_speed_in_s, ...
            %     config_vehicle.maximum_speed, config_vehicle.corner_speed_acc_cycle);
            dc = create_acc_cycle(10, 160, 100);
        case 'driving'
            dc = load("././data/DrivingCycles/"+config_vehicle.driving_cycle+".mat").dc;
        otherwise
            dc = load("././data/DrivingCycles/"+cycle_name+".mat").dc;
    end
    SharedConfig = evalin('base', 'SharedConfig');
    set_param(SharedConfig,'StopTime', num2str(length(dc.time)));
    disp("Simulation Stop at " + get_param(SharedConfig,'StopTime') + "s")
    assignin('base', "SharedConfig", SharedConfig);
    assignin('base', 'dc', dc);
end
