function Mode = updateInverterInformation(config_vehicle, config_motor)


     if config_motor.type1 == "ASM" && config_vehicle.inverter_type == "IGBT"
        Mode = 1;
        warning('off', 'Simulink:Commands:ParamUnknown');
    elseif config_motor.type1 == "ASM" && config_vehicle.inverter_type == "MOSFET"
        Mode = 2;
        warning('off', 'Simulink:Commands:ParamUnknown');
    elseif config_motor.type1 == "PMSM" && config_vehicle.inverter_type == "IGBT"
        Mode = 3;
        warning('off', 'Simulink:Commands:ParamUnknown');
    elseif config_motor.type1 == "PMSM" && config_vehicle.inverter_type == "MOSFET"
        Mode = 4;
        warning('off', 'Simulink:Commands:ParamUnknown');
     end
end