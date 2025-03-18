function [battery_GD, config_battery] = InitBattery(config_vehicle)
    % INITBATTERY - Initializes the battery configuration based on vehicle settings.
    %
    % Inputs:
    %   config_vehicle - Structure containing vehicle configuration details.
    %
    % Outputs:
    %   battery_GD     - Preallocated battery configuration structure.
    %   config_battery - Battery configuration loaded from JSON.
    
    % Load standard battery configuration from JSON file
    config_battery = read_json_config("config/conf_battery.json");
    
    % Assign vehicle configuration to the battery configuration
    config_battery.vehicle = config_vehicle;

    % Set default inverter type for the vehicle
    config_vehicle.inverter_type = 'MOSFET'; 

    % Additional battery properties
    config_battery.libraryname = "GD_Batterymodel_NSGAII_1"; % Name of created battery library
    
    % Dimension scaling factors for battery packaging
    config_battery.dim_x_factor = 0.75; 
    config_battery.dim_y_factor = 0.85; 
    config_battery.dim_z_factor = 0.2;

    % Preallocate battery configuration structure
    battery_GD = preallocate_configs_1_mod_all;
    
    % Initialize module and system information structures
    battery_GD.ModInfo = struct(); 
    battery_GD.SysInfo = struct(); 
end 

%% helper functions
function config_x = read_json_config(file_path)
    % Read json formatted table of config data and import it to Matlab workspace
    fid = fopen(file_path, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    config_x = jsondecode(str);
end