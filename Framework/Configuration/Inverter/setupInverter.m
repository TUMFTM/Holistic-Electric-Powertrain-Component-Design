function config_Inverter = setupInverter()

    config_Inverter = read_json_config("config/conf_Inverter.json");

    %Mode_IM_IGBT = Simulink.Variant('Mode == 1');
    %Mode_IM_MOSFET = Simulink.Variant('Mode == 2');
    %Mode_PSM_IGBT = Simulink.Variant('Mode == 3');
    %Mode_PSM_MOSFET = Simulink.Variant('Mode == 4');
   
    assignin('base', 'Mode_IM_IGBT', Simulink.Variant('Mode == 1'));
    assignin('base', 'Mode_IM_MOSFET', Simulink.Variant('Mode == 2'));
    assignin('base', 'Mode_PSM_IGBT', Simulink.Variant('Mode == 3'));
    assignin('base', 'Mode_PSM_MOSFET', Simulink.Variant('Mode == 4'));
end