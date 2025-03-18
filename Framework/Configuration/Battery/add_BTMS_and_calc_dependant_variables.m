function configs_6_BTMS_passed = add_BTMS_and_calc_dependant_variables(configs_4_sys_passed, user)
    configs_5_BTMS_all = preallocate_configs_5_BTMS;

    for ii = 1:1:size(configs_4_sys_passed, 2)
    
        config_bat = configs_4_sys_passed(ii);
        
        % Determine BTMS architecture
    
        config_bat.BTMSConfig = main_create_BTMS_config(config_bat.BTMSPara, config_bat.SysPara, config_bat.SysInfo, config_bat.ModInfo);
        config_bat.BTMS_ID = config_bat.BTMSPara.name;
        
        % Calculate thermal system properties depending on module and BTMS config
        
    %config.SysPara.thermal = calc_thermal_system_properties(config.BTMSConfig, config.SysPara.thermal);
        
        configs_5_BTMS_all = append_configs(configs_5_BTMS_all, config_bat);
    
    end
    configs_6_BTMS_passed = preallocate_configs_6_BTMS; % Preallocating the cell-array with all configuations that passed the module tests
    configs_6_BTMS_failed = preallocate_configs_6_BTMS; % Preallocating the cell-array with all configuations that failed the module tests
    
    for ii = 1:1:size(configs_5_BTMS_all, 2)
        
        config_bat = configs_5_BTMS_all(ii);
    
        % Exclude Configs that violate the max. dimensions (with consideration of BTMS)
    
        [config_bat, passed_BTMS] = test_dimensions_BTMS(config_bat);
    
        % Exclude configs that did not pass the tests from further consideration,
        % pass on working configs to the next steps.
    
        if check_for_failed_tests(passed_BTMS) && ~user % Test if any test has failed
            configs_6_BTMS_failed = append_configs(configs_6_BTMS_failed, config_bat, passed_BTMS, 'fail');
    
        else    % Those configurations have passed
            configs_6_BTMS_passed = append_configs(configs_6_BTMS_passed, config_bat, passed_BTMS, 'pass');
        end
    end

end