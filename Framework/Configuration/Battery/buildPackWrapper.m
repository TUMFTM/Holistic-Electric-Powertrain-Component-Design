function battery_GD = buildPackWrapper(battery_GD)
    % buildPackWrapper
    % Recalculates pack configurations and updates battery_GD.

    % Step 1: Recalculate system-level configurations
    battery_GD = parallelCellDistributionSys(battery_GD);

    % Step 2: Calculate system-dependent variables
    battery_GD = calc_sys_dependent_variables(battery_GD, true);

    % Step 3: Add BTMS and calculate dependent variables
    %battery_GD = add_BTMS_and_calc_dependant_variables(battery_GD, true);

    if numel(battery_GD) > 1
        warning('Multiple systems detected in battery_GD. Displaying information for the first system only.');
        selectedSystem = battery_GD(1); % Take the first system
    else
        selectedSystem = battery_GD; % Single system case
    end

    % Ensure SysInfo field exists in the selected system
    if isfield(selectedSystem, 'SysInfo')
        % Display the recalculated system parameters for the selected system
        fprintf('System Voltage: %.2f V\n', selectedSystem.SysInfo.U_nom_sys);
        fprintf('System Current: %.2f A\n', selectedSystem.SysInfo.I_max_sys);
        fprintf('System Capacity: %.2f Ah\n', selectedSystem.SysInfo.C_sys);
    end    
end