classdef init_simulation
    methods(Static)

        function batteryInitFlag = init_battery(bat_size_factor)
            %% Battery
            if nargin < 1 || isempty(bat_size_factor)
                bat_size_factor = 0;
            end
            config_vehicle = evalin('base', 'config_vehicle');
            config_battery = evalin('base', 'config_battery');
            if config_battery.E_sys_min + bat_size_factor > 10
                config_battery.E_sys_min = config_battery.E_sys_min + bat_size_factor;
            else
                config_battery.E_sys_min = 10;
            end

            % load pathes
            pathes.bat = jsondecode(fileread("input_and_parameters/00_programm_pathes/sim_BTMS_pathes.json"));
            pathes.bat = pathbuilder(pathes.bat, "models/battery/");

            battery = struct();
            battery = sim_BTMS_interface(battery,config_battery, config_vehicle); % call interface
            assignin('base','battery',battery)

            % Load battery config - make path more modulare with existing
            % path variable
            
            h.A = load(pathes.bat.Batterymodel+"/"+config_battery.libraryname+"/"+config_battery.libraryname+"_param_dch.mat");
            h.f = (fieldnames(h.A));
            clear(h.f{:})
            for i = 1:(size(h.f,1))
                assignin('base',string(h.f(i)),h.A.(string((h.f(i)))))
            end
            clear direction

            disp("init_battery: Successfully designed and initialized the battery module")
            batteryInitFlag = 1;
        end
          
    end
end