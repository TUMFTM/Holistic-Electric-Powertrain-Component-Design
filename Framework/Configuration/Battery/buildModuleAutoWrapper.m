function battery_GD = buildModuleAutoWrapper(battery_GD, input_configs)

    % Step 2: Recalculate configurations for module connections
    configs_1_mod_all = automaticModuleConnectionsCalculation(battery_GD, input_configs);

    % Step 3: Calculate dependent variables and update battery_GD
    battery_GD = calc_mod_dependent_variables(configs_1_mod_all, false);
    
    
    % Display the recalculated module parameters
    if numel(battery_GD) > 1
         %fprintf('Multiple modules detected in battery_GD. Displaying information for the first module only.');
    %    selectedModule = battery_GD(1); % Take the first module
        fprintf('Several modules have been successfully saved.')
    else
    %    selectedModule = battery_GD; % Single module case
         fprintf('A module have been succseffully saved.')
    end

    %% Display the recalculated module parameters for the selected module
    %fprintf('Module Voltage: %.2f V\n', selectedModule.ModInfo.U_nom_mod);
    %fprintf('Module Capacity: %.2f Ah\n', selectedModule.ModInfo.C_mod);
    %fprintf('Module Current: %.2f A\n', selectedModule.ModInfo.I_max_mod);
end