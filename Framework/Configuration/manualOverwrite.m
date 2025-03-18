function config_battery = manualOverwrite(config_battery, battery_GD)
    % Manually overwrite fields in config_battery with corresponding fields from battery_GD

    % Overwrite num_higher_p_mod
    if isfield(battery_GD, 'SysSpec') && isfield(battery_GD.SysSpec, 'num_higher_p_mod')
        config_battery.num_higher_p_mod = battery_GD.SysSpec.num_higher_p_mod;
    end

    % Overwrite num_higher_p_sys
    if isfield(battery_GD, 'SysSpec') && isfield(battery_GD.SysSpec, 'num_higher_p_sys')
        config_battery.num_higher_p_sys = battery_GD.SysSpec.num_higher_p_sys;
    end

    % Overwrite cellname (cell_id in battery_GD)
    if isfield(battery_GD, 'cell_ID')
        config_battery.cellname = battery_GD.cell_ID;
    end

    % Overwrite parameters from U_mod_nom to dim_z_sys_max
    if isfield(battery_GD, 'SysSpec')
        fieldsToOverwrite = {'U_mod_nom', 'U_mod_max', 'U_mod_min', ...
                             'U_sys_nom', 'U_sys_max', 'U_sys_min', ...
                             'C_sys_min', 'E_sys_min', 'dim_x_sys_max', ...
                             'dim_y_sys_max', 'dim_z_sys_max'};
                         
        for i = 1:length(fieldsToOverwrite)
            field = fieldsToOverwrite{i};
            if isfield(battery_GD.SysSpec, field)
                config_battery.(field) = battery_GD.SysSpec.(field);
            end
        end
    end

    % Overwrite Cell_capacity (C_A in battery_GD.BatPara.electrical)
    if isfield(battery_GD, 'BatPara') && isfield(battery_GD.BatPara, 'electrical') && isfield(battery_GD.BatPara.electrical, 'C_A')
        config_battery.Cell_capacity = battery_GD.BatPara.electrical.C_A;
    end
end