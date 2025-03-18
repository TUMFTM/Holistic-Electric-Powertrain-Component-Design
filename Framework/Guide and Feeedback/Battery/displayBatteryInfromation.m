function pack_info_table = displayBatteryInfromation(battery_GD)
    % Extract pack information from the structure
    capacity = battery_GD.SysInfo.C_sys;        % Pack capacity (Ah)
    energy = battery_GD.SysInfo.E_sys;          % Pack energy (kWh)
    nominal_voltage = battery_GD.SysInfo.U_nom_sys; % Nominal voltage (V)
    max_current = battery_GD.SysInfo.I_max_sys; % Maximum current (A)
    weight = battery_GD.SysInfo.mass_sys;       % Pack weight (kg)
    dim_x = battery_GD.SysInfo.dim_x_sys;       % Pack dimension (x) (m)
    dim_y = battery_GD.SysInfo.dim_y_sys;       % Pack dimension (y) (m)
    dim_z = battery_GD.SysInfo.dim_z_sys;       % Pack dimension (z) (m)

    % Convert dimensions from meters to millimeters for better readability
    dim_x_mm = dim_x * 1000; 
    dim_y_mm = dim_y * 1000; 
    dim_z_mm = dim_z * 1000;

    % Create a table to display the information
    pack_info_table = table( ...
        {'Capacity (Ah)'; 'Energy (kWh)'; 'Nominal Voltage (V)'; ...
         'Maximum Current (A)'; 'Dimension Length (mm)'; ...
         'Dimension Width (mm)'; 'Dimension Height (mm)'; 'Weight (kg)'}, ...
        [capacity; energy; nominal_voltage; max_current; ...
         dim_x_mm; dim_y_mm; dim_z_mm; weight], ...
        'VariableNames', {'Parameter', 'Value'} ...
    );

    % Display the table inline in the Live Script
    %disp('Battery Pack Information:');
    %disp(pack_info_table);
end