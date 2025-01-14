function [config, passed] = sys_test_current(config)

% This function tests the energy in kWh of a pack and compares it to the
% specifications made in 'system_specifications'

%% Calculate some module info and module energy
config.SysInfo.I_max_sys = config.BatPara.electrical.Crate_max * config.SysInfo.C_sys;          % Maximum module voltage in V

%% Test system against criteria

passed = struct('current_sys',false);

if config.SysInfo.I_max_sys >= config.SysSpec.I_sys_max
    passed.current_sys = true;
end
