function configs_3_sys_all = automaticPackConnectionsCalculation(battery_GD)            
             
    configs_3_sys_all = preallocate_configs_3_sys(); % Preallocate configurations array

    % Iterate through configurations in battery_GD
    for ii = 1:size(battery_GD, 2)
        config_bat = battery_GD(ii);

        % Calculate serial and parallel configurations
        s_sys = ceil(config_bat.SysSpec.U_sys_nom / config_bat.ModInfo.U_nom_mod);
        config_bat.SysPara.s_sys = s_sys;

        p_sys_raw = max(ceil(config_bat.SysSpec.I_sys_max / config_bat.ModInfo.I_max_mod), ...
                        ceil(config_bat.SysSpec.C_sys_min / config_bat.ModInfo.C_mod));
        config_bat.SysPara.p_sys = p_sys_raw;
      
        % Determine spatial arrangements
        epe_sys = epe_distribution(p_sys_raw, config_bat.SysSpec.num_higher_p_sys, ...
                                    config_bat.SysSpec.sys_min_e, config_bat.SysSpec.sys_max_e);
        for jj = 1:size(epe_sys, 2)
            p_sys = epe_sys(jj).p;
            for kk = 1:size(epe_sys(jj).pe, 2)
                % Populate system parameters
                config_bat.sys_ID = get_config_ID(); % Assign unique ID
                config_bat.SysInfo.num_mods_sys = s_sys * p_sys;
                config_bat.SysInfo.num_serial_mods_sys = s_sys;
                config_bat.SysInfo.num_parallel_mods_sys = p_sys;
                config_bat.SysInfo.num_layers_sys = epe_sys(jj).e(kk);
                config_bat.SysInfo.num_parallel_mods_per_layer_sys = epe_sys(jj).pe(kk);


                config_bat.SysPara.p_sys = p_sys * config_bat.SysPara.p_mod;
                config_bat.SysPara.s_sys = s_sys * config_bat.SysPara.s_mod;
                config_bat.SysPara.pe_sys = epe_sys(jj).pe(kk) * config_bat.SysPara.pe_mod;
                config_bat.SysPara.e_sys = epe_sys(jj).e(kk) * config_bat.SysPara.e_mod;

                %fieldsToRemove = {'BTMS_ID', 'BTMSConfig', 'Tests_BTMS'};

                % Loop through and remove fields if they exist
                %for i = 1:numel(fieldsToRemove)
                %    if isfield(config_bat, fieldsToRemove{i})
                %        config_bat = rmfield(config_bat, fieldsToRemove{i});
                %    end
                %end

                %if isfield(config_bat, 'Tests_sys')
                %    config_bat = rmfield(config_bat, 'Tests_sys');
                %end

                configs_3_sys_all = append_configs(configs_3_sys_all, config_bat); % Append configuration
            end
        end
    end
end    