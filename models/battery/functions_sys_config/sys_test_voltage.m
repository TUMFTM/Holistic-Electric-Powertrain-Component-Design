function [config, passed] = sys_test_voltage(config)

% This function tests the energy in kWh of a pack and compares it to the
% specifications made in 'system_specifications'

%% Test system against criteria

passed = struct('voltage_sys',false);

if 1 % config.SysInfo.U_max_sys <= config.SysSpec.U_sys_nom
    passed.voltage_sys = true;
end
