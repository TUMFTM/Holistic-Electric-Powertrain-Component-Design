function import_dc_cycle(cycle_name)
    config_vehicle = evalin('base', 'config_vehicle');
    switch cycle_name
        case 'custom'
            dc = create_acc_cycle(config_vehicle.zero_to_corner_speed_in_s, ...
                config_vehicle.maximum_speed, config_vehicle.corner_speed_acc_cycle);
        case 'driving'
            dc = load("././data/DrivingCycles/"+config_vehicle.driving_cycle+".mat").dc;
        otherwise
            dc = load("././data/DrivingCycles/"+cycle_name+".mat").dc;
    end

    disp("Simulation Stop at " + num2str(dc.time(end)) + "s")
    assignin('base', 'dc', dc);
end
