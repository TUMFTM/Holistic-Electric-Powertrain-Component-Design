function module_info_table = displayModuleInformation(battery_GD)
    
    if numel(battery_GD) > 1
        fprintf('Multiple modules detected in battery_GD. Displaying information for the first module only.');
        battery_GD = battery_GD(1); % Take the first module
    end
    % Extract module information from the structure
    capacity = battery_GD.ModInfo.C_mod;        % Module capacity
    nominal_voltage = battery_GD.ModInfo.U_nom_mod; % Nominal voltage
    max_current = battery_GD.ModInfo.I_max_mod; % Maximum current
    weight = battery_GD.ModInfo.mass_mod;       % Module weight
    dim_x = battery_GD.ModInfo.dim_x_mod;     % Module dimension (x)
    dim_y = battery_GD.ModInfo.dim_y_mod;       % Module dimension (y)
    dim_z = battery_GD.ModInfo.dim_z_mod;       % Module dimension (z)
   

    % Create a table to display the information
     module_info_table = table( ...
        {'Capacity (Ah)'; 'Nominal Voltage (V)'; 'Maximum Current (A)'; ...
         'Dimension Length (m)'; 'Dimension Width (m)'; 'Dimension Height (m)'; 'Weight (kg)'}, ...
        [capacity; nominal_voltage; max_current; dim_x; dim_y; dim_z; weight], ...
        'VariableNames', {'Parameter', 'Value'} ...
    );



    % Display the table inline in the Live Script
    %disp('Battery Module Information:');
    %disp(module_info_table);
end