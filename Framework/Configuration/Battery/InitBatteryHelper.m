function [battery_GD, SysSpec, input_configs] = InitBatteryHelper(BatPara, config_battery, battery_GD, config_vehicle)
    % INITBATTERYHELPER - Initializes battery system configurations based on cell type.
    %
    % Inputs:
    %   BatPara        - Struct containing battery parameters.
    %   config_battery - Struct containing battery configuration.
    %   battery_GD     - Struct for battery global data.
    %   config_vehicle - Struct containing vehicle configuration.
    %
    % Outputs:
    %   battery_GD     - Updated battery global data.
    %   SysSpec        - System specifications.
    %   input_configs  - Configuration inputs for the system.

    % ------------------ DEVELOPMENT CONTROL ------------------
    % Skip testing stages if needed (set testing_skipper > 0 during development)
    testing_skipper = 0;

    % ------------------ STEP 1: CONFIGURE CELL TYPE ------------------
    if testing_skipper < 2
        try
            switch lower(BatPara.cell_type)
                case {"cyl", "cylindrical"}
                    % Cylindrical cell configuration
                    BatPara.cell_type = "Cylindrical";
                    input_configs = {
                        config_battery.cellname, 'system_para_BTMS_sim', 'GD_liquid_inside_Cyl'; ...
                    };
                    run GD_system_specifcations_Cyl_inside;

                case "pouch"
                    % Pouch cell configuration
                    input_configs = {
                        config_battery.cellname, 'system_para_BTMS_sim', 'GD_liquid_inside_Pouch'; ...
                    };
                    run GD_system_specifcations_Pouch_inside;

                case "prismatic"
                    % Prismatic cell configuration
                    input_configs = {
                        config_battery.cellname, 'system_para_BTMS_sim', 'GD_liquid_inside_Prismatic'; ...
                    };
                    run GD_system_specifcations_Prismatic_inside;

                otherwise
                    error("Unsupported cell type: %s. Please check the configuration file.", BatPara.cell_type);
            end

            % Enforce system specifications based on selected configuration
            SysSpec = config_enforcer(SysSpec, config_battery); %#ok<NODEF>

        catch ME
            % Error handling for cell configuration issues
            error("Error in cell configuration: %s", ME.message);
        end
    end

    % ------------------ STEP 2: ADJUST SYSTEM DIMENSIONS ------------------
    try
        % Scale battery system dimensions based on vehicle size
        SysSpec.dim_x_sys_max = config_battery.dim_x_factor * config_vehicle.vehicle_dim_x_mm;
        SysSpec.dim_y_sys_max = config_battery.dim_y_factor * config_vehicle.vehicle_dim_y_mm;
        SysSpec.dim_z_sys_max = config_battery.dim_z_factor * config_vehicle.vehicle_dim_z_mm;

        % Assign updated system specifications to battery_GD
        battery_GD.SysSpec = SysSpec;
        
    catch ME
        % Error handling for system dimension adjustments
        error("Error adjusting system dimensions: %s", ME.message);
    end
end