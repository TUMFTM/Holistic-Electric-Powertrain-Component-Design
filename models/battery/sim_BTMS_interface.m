function [battery_GD] = sim_BTMS_interface (battery_GD,config_battery, config_vehicle)

%% sim_BTMS_interface v0.15
% This interface for sim_BTMS was created as part of Global Drive 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                           Versions                                    %
% v0.1 - first setup of System Parametrisation                          %
% v0.2 - automatetd choice of cooling solution based on battery cells   %
% v0.3 - adaption towards GD-framework                                  %
% v0.4 - added config_enforcer                                          %
% v0.5 - added pathbuilder                                              %
% v0.6 - changed sim_BTMS_interface into function                       %
% v0.7 - running version with complete static calculation               %
% v0.8 - added simscape battery - build battery                         %
% v0.9 - implemented Reiter Simulation                                  %
% v0.10 - simscape battery with charge dynamics                         %
% v0.11 - adaption towards msimulationframework                         %
% v0.12 - Designs battery based on max allowed nominal voltage          %
% v0.13 - added ManufacturerParameterized Cells                         %
% v0.14 - added LossProbeAdder                                          %
% v0.14 - integration of build cells                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        known Problems                                 %
% - Skipper broken                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Testparameters

% Skipper for faster testing
testing_skipper = 0;

if testing_skipper > 0
    fprintf("!!!\n!!!\n\nSkipping sim_BTMS_1 directly to step %d\n\n!!!\n!!!\n", testing_skipper)
    configs_6_BTMS_passed = "configs_6_BTMS_passed.mat";
    load(configs_6_BTMS_passed)
    battery_GD.config = configs_elector(battery_GD, configs_6_BTMS_passed);
end    

%% 1 Initialation
import simscape.battery.builder.*;

if testing_skipper < 1

%Load GD battery config 
config_battery.vehicle = config_vehicle;
end

% load pathes
pathes.bat = jsondecode(fileread("input_and_parameters/00_programm_pathes/sim_BTMS_pathes.json"));
pathes.bat = pathbuilder(pathes.bat, "models/battery/");

BatPara = struct();

%% Debugging/Testing Stuff

%Visualization of the build battery
visualization.Cell = 0;
visualization.Module_P = 0;
visualization.Module_S = 0;
visualization.Module_A = 0;
visualization.Pack = 0;

%% 2 Interpret Configuration

if testing_skipper < 2
    if config_battery.manufacturer_parameterization == 0
        % For Cells we provided the data for
        BatPara = get_BatPara(BatPara, config_battery);

    elseif ~strcmp(config_battery.manufacturer, 'None')
        % For cells included in the Simscape Package
        
        batteryCell = Cell;
        
        % Check if chosen Cell exists
        h.manufacturers_available = batteryCell.AvailableParameterizationManufacturers;
        if ismember(config_battery.manufacturer, h.manufacturers_available)
            batteryCell.ParameterizationManufacturer = config_battery.manufacturer;
            h.parameterization_available = batteryCell.AvailableParameterizationPartNumbers;
            if ismember(config_battery.cellname, h.parameterization_available)
                batteryCell.ParameterizationPartNumber = config_battery.cellname;
                
                [BatPara] = CellParaBuilder(BatPara, config_battery);
                % Deactivate Charge Dynamics for PreParameterized Cells
                config_battery.simulation_charge_dynamics = 'rc0';
                disp("Turned off Charge Dynamics (->RC0) as a PreParameterized batterycell was chosen");
            else
                error("Not allowed cellinformation! Unknown cellname")
            end
        else
            error("Not allowed cellinformation! Unknown manufacturer")
        end
    else
        error("Not allowed cellinformation! Set manufacturer_parameterization to 0 or choose a valid manufacturer")
    end
    
    % Choosing configuration based on cell_type                
        if ismember(BatPara.cell_type, ["Cyl", "Cylindrical"])
            BatPara.cell_type = "Cylindrical";
            input_configs = {
            config_battery.cellname,  'system_para_BTMS_sim', 'GD_liquid_inside_Cyl'; ...
            };
            run GD_system_specifcations_Cyl_inside
        
        elseif BatPara.cell_type == "Pouch"
        
            input_configs = {
            config_battery.cellname,  'system_para_BTMS_sim', 'GD_liquid_inside_Pouch'; ...
            };
            run GD_system_specifcations_Pouch_inside

        elseif BatPara.cell_type == "Prismatic"
        
            input_configs = {
            config_battery.cellname,  'system_para_BTMS_sim', 'GD_liquid_inside_Prismatic'; ...
            };
            run GD_system_specifcations_Prismatic_inside
        
        else
            error("unsupported celltype/n check configuration file")
            
        end
        
        % Enforce GD config values
        SysSpec = config_enforcer(SysSpec, config_battery); %#ok<NODEF>

end

%% 3 Create and choose battery system and BTMS concept

if testing_skipper < 3

% call sim_BTMS_1_setup to creat battery systems and BTMS concepts
% sim_BTMS_1_caller(BatPara, config, input_configs, pathes, SysSpec)
configs_6_BTMS_passed = main_sim_BTMS_1_system_setup(BatPara, config_battery, input_configs, pathes, SysSpec);
clearvars -except battery_GD pathes config_battery configs_6_BTMS_passed visualization testing_skipper

% Calculate cost
%%%% To be added by Gideon %%%%

% Cost in € - Placeholders
% for i = 1:size(configs_6_BTMS_passed,2)
%     configs_6_BTMS_passed(i).SysInfo.Cost = costdata_evaluator(configs_6_BTMS_passed(i));
%     %configs_6_BTMS_passed(i).SysInfo.Cost = costcode_bat(configs_6_BTMS_passed(i));
%     %battery_GD.SysInfo.Cost.Germany = 11.800; 
%     %battery_GD.SysInfo.Cost.USA = 11.000;
%     %battery_GD.SysInfo.Cost.China = 10.800;
% end

% chooses best fitting battery configuration
battery_GD = configs_elector(battery_GD, configs_6_BTMS_passed);
clearvars configs_6_BTMS_passed % clear not chossen configurations


% Cost in € !!!Caution, cost module only working with Windows!!!!
% battery_GD.SysInfo.Cost = costdata_evaluator(battery_GD);

% Calculate Sustainability Index
battery_GD = sustainability_estimator(battery_GD);

end

%% 4 Simscape Battery

if  config_battery.simulation_simscape == 1 && testing_skipper < 4
    battery_builder(battery_GD, config_battery, pathes, visualization);   % build battery pack
    battery_GD = battery_initializer(battery_GD, config_battery, pathes);    % initialize battery pack
    battery_PackReplacement = battery_PackReplacer; % Ensures usage of the newly created battery pack
    battery_LossProbe = battery_LossProbeAdder(battery_GD); % Add Probe to get power dissipation
    %run batterySimulation
end

%% 5 Simulation Reiter

if config_battery.simulation_reiter == 1 && testing_skipper < 4
    main_sim_BTMS_2_system_simulation(battery_GD, config_battery)
end

%% Help Functions

function [BatPara] = get_BatPara(BatPara, config)
    t_cell_capacity = 0;
    try
        run (config.cellname);
    catch
        load("models/battery/input_and_parameters/01_cell_data/"+ config.cellname + "/" + config.cellname+ ".mat")
    end
end

%function [] = sim_BTMS_1_caller(BatPara, config, input_configs, pathes, SysSpec)
%    run main_sim_BTMS_1_system_setup
%end

end
