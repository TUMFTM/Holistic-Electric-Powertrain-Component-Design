function configs_3_sys_all = parallelCellDistributionSys(configs_2_mod_passed)
    % parallelCellDistributionSys: Calculates and organizes the spatial 
    % arrangement of parallel connections in the battery system.
    %
    % Functionality:
    %   - Distributes modules in a 3D spatial arrangement (e * pe) for parallel connections.
    %   - Updates system information (`SysInfo`) and system parameters (`SysPara`) for each configuration.
    %   - Appends valid configurations to a cell array for further processing.
    %
    % Input:
    %   - config_bat: A structure containing the current battery configuration, 
    %                 including system parameters (`SysPara`) and specifications (`SysSpec`).
    %
    % Output:
    %   - config_bat: The updated configuration after processing parallel connections.
    
    % Preallocate the configurations array
    configs_3_sys_all = preallocate_configs_3_sys();
    
    
    for ii = 1:size(configs_2_mod_passed, 2)
        config_bat = configs_2_mod_passed(ii);
        % Generate spatial arrangements for parallel connections
        % `epe_sys` represents the possible combinations of e (layers) and pe (parallel modules per layer).
        epe_sys = epe_distribution(config_bat.SysInfo.num_parallel_mods_sys, ...
                                    config_bat.SysSpec.num_higher_p_sys, ...
                                    config_bat.SysSpec.sys_min_e, ...
                                    config_bat.SysSpec.sys_max_e);
        s_sys = config_bat.SysInfo.num_serial_mods_sys;
        % Iterate through all possible parallel connection configurations (p)
        for jj = 1:size(epe_sys, 2)
           
            p_sys = epe_sys(jj).p; % Number of parallel modules
           
            % Iterate through all combinations of layers (e) and modules per layer (pe)
            for kk = 1:size(epe_sys(jj).pe, 2)
                % Update system information with current configuration
                config_bat.SysInfo.num_mods_sys = s_sys * p_sys; % Total modules
                config_bat.SysInfo.num_layers_sys = epe_sys(jj).e(kk); % Number of layers
                config_bat.SysInfo.num_parallel_mods_sys = p_sys;
                config_bat.SysInfo.num_parallel_mods_per_layer_sys = epe_sys(jj).pe(kk); % Modules per layer
    
                % Update system parameters
                config_bat.SysPara.p_sys = p_sys * config_bat.SysPara.p_mod; % Total parallel modules
                config_bat.SysPara.s_sys = s_sys * config_bat.SysPara.s_mod; % Total serial modules
                config_bat.SysPara.pe_sys = epe_sys(jj).pe(kk) * config_bat.SysPara.pe_mod; % Parallel modules per layer
                config_bat.SysPara.e_sys = epe_sys(jj).e(kk) * config_bat.SysPara.e_mod; % Number of layers
    
                % Assign a unique configuration ID
                config_bat.sys_ID = get_config_ID();
    
                % Remove `Tests_sys` field to avoid errors in downstream processing
               %if isfield(config_bat, 'Tests_sys')
               %     config_bat = rmfield(config_bat, 'Tests_sys');
                %end
                %if isfield(config_bat, 'BTMS_ID')
                %    config_bat = rmfield(config_bat, 'BTMS_ID');
                %end
                %if isfield(config_bat, 'BTMSConfig')
                %    config_bat = rmfield(config_bat, 'BTMSConfig');
                %end
                %if isfield(config_bat, 'Tests_BTMS')
                %    config_bat = rmfield(config_bat, 'Tests_BTMS');
                %end
    
                
                % Append the current configuration to the list of all configurations
                configs_3_sys_all = append_configs(configs_3_sys_all, config_bat);
            end
        end
    end   
end