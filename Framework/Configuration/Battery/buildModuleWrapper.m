function battery_GD = buildModuleWrapper(battery_GD, input_configs)

    % Step 2: Recalculate configurations for a single module
    configs_1_mod_all = parallelCellDistribution(battery_GD, input_configs);

    % Step 3: Calculate dependent variables and update battery_GD
    try
        battery_GD = calc_mod_dependent_variables(configs_1_mod_all, true);

        % Validate the output
        if ~isstruct(battery_GD)
            error('InvalidOutput:NotStruct', 'The variable battery_GD is not a structure.');
        elseif isempty(battery_GD)
            error('InvalidOutput:EmptyStructure', 'The variable battery_GD is empty.');
        end

        % Check if `battery_GD` has multiple elements
        if numel(battery_GD) > 1
            warning('Multiple modules detected in battery_GD. Displaying information for the first module only.');
            selectedModule = battery_GD(1); % Take the first module
        else
            selectedModule = battery_GD; % Single module case
        end

        % Display the recalculated module parameters for the selected module
        fprintf('Module Voltage: %.2f V\n', selectedModule.ModInfo.U_nom_mod);
        fprintf('Module Capacity: %.2f Ah\n', selectedModule.ModInfo.C_mod);
        fprintf('Module Current: %.2f A\n', selectedModule.ModInfo.I_max_mod);

    catch ME
        % Handle and display errors gracefully
        fprintf(2, 'Error: %s\n', ME.message); % Print error message in red
        switch ME.identifier
            case 'InvalidOutput:NotStruct'
                fprintf(2, 'Please ensure calc_mod_dependent_variables returns a structure.\n');
            case 'InvalidOutput:EmptyStructure'
                fprintf(2, 'The variable battery_GD is empty. Please check the input configurations.\n');
            otherwise
                fprintf(2, 'An unexpected error occurred. Please debug the function.\n');
        end
    end
end