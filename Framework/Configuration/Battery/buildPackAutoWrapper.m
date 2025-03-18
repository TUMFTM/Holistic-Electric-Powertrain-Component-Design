function [battery_GD, f] = buildPackAutoWrapper(battery_GD, user)
    % buildPackAutoWrapper
    % Automatically calculates and updates pack-level configurations with error handling.

    try
        % Step 1: Recalculate pack-level configurations
        configs_3_sys_all = automaticPackConnectionsCalculation(battery_GD);

        % Step 2: Calculate system-dependent variables
        battery_GD = calc_sys_dependent_variables(configs_3_sys_all, user);

        % Validate the configurations after system-dependent calculations
        %validateBatteryConfigurations(battery_GD, 'pack-level parameters');

        % Step 3: Add BTMS and calculate dependent variables
        %battery_GD = add_BTMS_and_calc_dependant_variables(battery_GD, false);

        % Validate configurations again after BTMS addition
        %validateBatteryConfigurations(battery_GD, 'BTMS parameters');

        % Step 4: Display final battery configurations
        [battery_GD, f] = displayBatteryConfigurations(battery_GD);

        fprintf('System Voltage: %.2f V\n', battery_GD.SysInfo.U_nom_sys);
        fprintf('System Current: %.2f A\n', battery_GD.SysInfo.I_max_sys);
        fprintf('System Capacity: %.2f Ah\n', battery_GD.SysInfo.C_sys);

    catch ME
        % Handle and display errors gracefully
        fprintf(2, 'Error: %s\n', ME.message); % Print error message in red
        if strcmp(ME.identifier, 'BatteryValidation:InvalidConfig')
            fprintf(2, 'Please verify the input configurations or adjust constraints.\n');
        else
            fprintf(2, 'An unexpected error occurred. Please debug the function.\n');
        end
    end
end

function validateBatteryConfigurations(battery_GD, context)
    % validateBatteryConfigurations
    % Checks for invalid battery configurations and throws an error if any are found.

    if size(battery_GD, 1) == 1 && any(isnan(battery_GD(1).mod_ID(:)))
        error('BatteryValidation:InvalidConfig', ...
              'No configurations were found. No %s were added to battery_GD.', context);
    end
end