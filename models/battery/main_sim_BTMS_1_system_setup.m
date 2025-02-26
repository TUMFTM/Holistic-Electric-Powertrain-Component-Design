function [configs_6_BTMS_passed] = main_sim_BTMS_1_system_setup(BatPara, config, input_configs, pathes, SysSpec)
    arguments
        BatPara struct = struct();
        config struct = struct();
        input_configs cell = cell();
        pathes struct = struct();
        SysSpec struct = struct();
    end

%% Info

% Part 1 of the BTMS configuration: Creation of the battery systems and
% BTMS concepts

% For analysis and simulation of the systems refer to 'main_sim_BTMS_2_system_simulation'.

%% GD

% sim_BTMS_1_system_setup v1.1
% Adapted for use in Global Drive 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                           Versions                                    %
% v1.1 - changed into functions to increase usability                   %
% v1.2 - Now chooses packconfiguration based on U_max                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        known Problems                                 %
% - can't run on its own anymore                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialation

% Check for usage inside the Global Drive Project
if exist('config','var') == 0
%if config.battery.name == [] % maybe implement a better solution later
    clearvars
    close all
    clc

    % loading save pathes - for GD this is done inside the sim_BTMS_interface
    pathes.bat = jsondecode(fileread('input_and_parameters/00_programm_pathes/sim_BTMS_pathes.json'));

    %% User Input: Cells and BTMS-method combinations in cosideration
    


    % Specify LIB parameters and corresponding BTMS method in the following format:

    % configurations = {
%     'Cell_1', 'SysInfo_1', 'BTMS_1';
%     'Cell_2', 'SysInfo_2', 'BTMS_2';
%     'Cell_3', 'SysInfo_3', 'BTMS_3';
%     'Cell_4', 'SysInfo_4', 'BTMS_4';
%     'Cell_5', 'SysInfo_5', 'BTMS_5';
%     };

    % Every row corresponds to a variant that gets considered for the battery
    % system configuration. The input data should be saved in
    % 'input_and_parameters/01_cell_data', '/02_BTMS_configs' and '/03_system_data', respectively.


    % % Cylindric and Pouch
% 
% input_configs = {
%     'Cyl_18650',  'system_para_BTMS_sim', 'liquid_bottom'; ...
%     'Cyl_21700',  'system_para_BTMS_sim', 'liquid_bottom'; ...
%     'Pouch_74Ah', 'system_para_BTMS_sim', 'liquid_bottom'; ...
%     'Cyl_18650',  'system_para_BTMS_sim', 'liquid_bottom_xdir'; ...
%     'Cyl_21700',  'system_para_BTMS_sim', 'liquid_bottom_xdir'; ...
%     'Pouch_74Ah', 'system_para_BTMS_sim', 'liquid_bottom_xdir'; ...
%     'Cyl_18650',  'system_para_BTMS_sim', 'liquid_full_system'; ...
%     'Cyl_21700',  'system_para_BTMS_sim', 'liquid_full_system'; ...
%     'Pouch_74Ah', 'system_para_BTMS_sim', 'liquid_full_system'; ...
%     };
% 
% % Provide the requested system specification
% run input_and_parameters\04_system_specifications\system_specifcations_cyl_Pouch.m;


    % % Prismatic (PHEV2)
% 
% % Two parallel serial module interconnections
% 
% input_configs = {    
%     'Pris_PHEV2', 'system_para_BTMS_sim', 'liquid_bottom'; ...   
%     'Pris_PHEV2', 'system_para_BTMS_sim', 'liquid_bottom_xdir'; ...   
%     'Pris_PHEV2', 'system_para_BTMS_sim', 'liquid_full_system'; ...
%     };
% 
% % Provide the requested system specification
% run input_and_parameters\04_system_specifications\system_specifcations_Pris_PHEV2.m;


    % % Prismatic (BEV2)
% 
% % Three parallel serial module interconnections 
% 
% input_configs = {    
%     'Pris_BEV2',  'system_para_BTMS_sim', 'liquid_bottom'; ...
%     'Pris_BEV2',  'system_para_BTMS_sim', 'liquid_bottom_xdir'; ...
%     'Pris_BEV2',  'system_para_BTMS_sim', 'liquid_full_system'; ...
%     };
% 
% % Provide the requested system specification
% run input_and_parameters\04_system_specifications\system_specifcations_Pris_BEV2.m;


    % % Cylindric (inside cooling)
    % 
    % % Module simulation with cooing inside the module
    % 
    % input_configs = {
    %     'Cyl_18650',  'system_para_BTMS_sim', 'liquid_inside_18650'; ...
    %     'Cyl_21700',  'system_para_BTMS_sim', 'liquid_inside_21700'; ...
    %     };
    % 
    % % Provide the requested system specification
    % run input_and_parameters\04_system_specifications\system_specifcations_cyl_inside.m;


    % % Pouch (inside cooling)
    %
    % Module simulation with cooing inside the module
    %
    %input_configs = {
    %    'Pouch_74Ah',  'system_para_BTMS_sim', 'liquid_inside_Pouch'; ...
    %    };

    % Provide the requested system specification
    %run input_and_parameters\04_system_specifications\system_specifcations_Pouch_inside.m;

end

%% Step 1: Creating the modules

% This section iterates through each configuration and creates the basic
% information about the modules. In the next step the modules are tested
% against the relevant criteria.

% Every possible configuration is stored in an array with all auxilary data
% needed for simulation. This is course creates somewhat of an overhead,
% but makes it easier to just pick up a configuration and start simulating.

clear append_configs_step_1                         % Clear persistent variable in function
configs_1_mod_all = preallocate_configs_1_mod_all;  % Preallocating the cell-array with all possible configurations


% Iterate through the configurations

for ii = 1:size(input_configs, 1)
    
    % Load cell data. 
    
    % This works, because we only provide the name of the dataset in
    % input_configs, so the parameters are loaded with a script.
    
    if config.manufacturer_parameterization == 0
        try
        run(input_configs{ii,1});
        catch
        load("models/battery/input_and_parameters/01_cell_data/"+ input_configs{ii,1} + "/" + input_configs{ii,1} + ".mat")
        end
    else

    end
    % Determine max size serial connection on module level
    
    %{
    % Take minimum of either specified maximum number of serially connected
    % cells or the maximum inplied by nominal module and nom. cell voltage
    
    s_mod_raw = min(SysSpec.s_mod_max, ceil(SysSpec.U_mod_nom / BatPara.electrical.U_nom));
    %}

    % Take minimum of either specified maximum number of serially connected
    % cells or the maximum inplied by nominal module and max cell voltage
    s_mod_raw = min(SysSpec.s_mod_max, floor(SysSpec.U_mod_nom / BatPara.electrical.U_max));
    
    % We only want to consider an even number of cells as an option -->
    % Round to the next even integer.
    
    s_mod = 2*floor(s_mod_raw/2);
    % s_mod = 2*ceil(s_mod_raw/2);
    
    % Determine min size of parallel connection
    
    % 
    % This is either defined by the max. individual cell current and the 
    % requested fast-charging capability (including 10% safety) or the module capacity.
    
    p_min_mod = max(ceil(SysSpec.I_mod_max * (1+SysSpec.sf_parallel) / BatPara.electrical.I_max), ceil(SysSpec.C_mod_min / BatPara.electrical.C_A));
   
    
    % Spatial arrangement of the parallel connection
    
    % We not only can arrange the cells in a simple one dimensional 's*p'
    % grid, but can arrange the cells in all three dimensions. Therefore in
    % this step we create an 'e*pe' connection from the p cells in
    % parallel so that p = e*pe. Those are arranged next to each other so
    % the total number of cells in the module is n = e*pe*s.
    
    epe_mod = epe_distribution(p_min_mod, SysSpec.num_higher_p_mod, SysSpec.mod_min_e, SysSpec.mod_max_e);   
   
    % Load BTMS data
        
    % This works, because we only provide the name of the dataset in
    % input_configs, so the parameters are loaded with a script.
    
    run(input_configs{ii,3});
    
    
    % More parallel connections may be considered for a given configuration
    % --> Interrate through all p's
    
    for jj = 1:1:size(epe_mod,2)   
        
        % To comply to https://github.com/TUMFTM/sim_battery_system as much
        % as possible we resuse the specifications of the serial and
        % parallel connection in 'system_data'
        
        % Load system data
                
        % Create 'p' and 's' fields to comply to naming convention from
        % 'sim_battery_system'.
        
        p = epe_mod(jj).p;  
        s = s_mod;
        
        run(input_configs{ii,2});
        
        % Delete 'p' and 's', we name them 'p_mod' and 's_mod' from now on,
        % to help with distinction from the system level.
        
        SysPara = rmfield(SysPara, {'p','s'});
        
        SysPara.p_mod = p;
        SysPara.s_mod = s;
        
        % Calc additional thermal cell parameters
        
        [SysPara.thermal.transfer.A_x, SysPara.thermal.transfer.A_y, SysPara.thermal.transfer.A_z, ...
            BatPara.thermal.A] = calc_cell_surface(BatPara);
        
        % For every p several combinations of e-pe may be found. --> We
        % iterate through them as well
        
        for kk = 1:1:size(epe_mod(jj).pe, 2)
            
            SysPara.pe_mod = epe_mod(jj).pe(kk);
            SysPara.e_mod = epe_mod(jj).e(kk);
          
            configs_1_mod_all = append_configs_step_1(BatPara, BTMSPara, SysPara, SysSpec, configs_1_mod_all); % Save everything in a cell-array 
        end
    end
end 


% Check if feasible configurations have been found 

check_feasible_configs(configs_1_mod_all,0)


% Save resuls of this step

%save("BTMS_simulation_results\configs_1_mod_all", 'configs_1_mod_all')  % Save all possible configurations of this run in a MAT-File
save(pathes.bat.save_1, 'configs_1_mod_all')  % Save all possible configurations of this run in a MAT-File

clearvars -except config* pathes* % Clear everything instead the array with the configs.

%% Step 2: Basic testing on module-level

% In this step we take the configurations from the last sections and test
% them for some basic critera. Here we only perform basic tests to decrease
% the solutions space as much as possible before we start with the
% computationally heavier stuff.

clear append_configs        % Clear persistent variable in function
load('configs_1_mod_all');  % Load the configurations from last step if not already in workspace

configs_2_mod_passed = preallocate_configs_2_mod; % Preallocating the cell-array with all configuations that passed the module tests
configs_2_mod_failed = preallocate_configs_2_mod; % Preallocating the cell-array with all configuations that failed the module tests


% Iterate through the configurations

for ii = 1:1:size(configs_1_mod_all, 2)
    
    config_bat = configs_1_mod_all(ii);

    % Exclude Configs that violate the max. dimensions (without consideration of BTMS)

    [config_bat, passed_dim] = mod_test_dimensions_no_BTMS(config_bat);


    % Exclude Configs that violate the energy content criteria

    [config_bat, passed_energy] = mod_test_energy(config_bat);


    % Create joint structure of criteria

    passed_mod = join_passed_structs(passed_dim, passed_energy);


    % Exclude configs that did not pass the tests from further consideration,
    % pass on working configs to the next steps.

    if check_for_failed_tests(passed_mod) % Test if any test has failed
        configs_2_mod_failed = append_configs(configs_2_mod_failed, config_bat, passed_mod, 'fail');

    else    % Those configurations have passed
        configs_2_mod_passed = append_configs(configs_2_mod_passed, config_bat, passed_mod, 'pass');
    end
end

% Check if feasible configurations have been found 

check_feasible_configs(configs_2_mod_passed,configs_2_mod_failed);


% Save results of this step

save(pathes.bat.save_2_passed, 'configs_2_mod_passed')  % Save all possible configurations of this run in a MAT-File
save(pathes.bat.save_2_failed, 'configs_2_mod_failed')  % Save all failed configurations of this run in a MAT-File
%save('BTMS_simulation_results\configs_2_mod_passed', 'configs_2_mod_passed')  % Save all possible configurations of this run in a MAT-File
%save('BTMS_simulation_results\configs_2_mod_failed', 'configs_2_mod_failed')  % Save all failed configurations of this run in a MAT-File

clearvars -except config* pathes* % Clear everything instead the array with the configs.


%% Step 3: Building a battery pack from the modules

% Now we combine the module configuration that passed the prior step to a
% battery pack. Basically the strategy is the same as in step 1: First we
% determine the number of serial and parallel module connections, then the
% arrange the modules connected in parallel spatially. 

% Note we don't have the BTMS heat-sinks included in the modules yet. Doing
% so would massively increase the number of possible configurations, so we
% first creates the possible battery packs and then check with which of our 
% modules we can fulfill our specified battery pack dimensions and energy
% criteria. What remains will be provided with a BTMS.

clear append_configs            % Clear persistent variable in function
clear get_config_ID             % Clear persistent variable in function
load('configs_2_mod_passed');   % Load the feasible configurations from last step if not already in workspace

configs_3_sys_all = preallocate_configs_3_sys; % Preallocating the cell-array with all configuations that passed the module tests


% Iterate through the configurations

for ii = 1:size(configs_2_mod_passed, 2)
    
    config_bat = configs_2_mod_passed(ii);

    % Determine number of modules connected in serial inside the battery system
    
    % U_max can be above config level - U_nominal at config level
    % This depends on the nominal module voltage and the specified nominal
    % system voltage.
     
    s_sys = ceil(config_bat.SysSpec.U_sys_nom / config_bat.ModInfo.U_nom_mod);
    
    %{ % U_max stays below config level
    % This depends on the maxiaml module voltage and the specified nominal
    % system voltage.
     
    % s_sys = floor(config_bat.SysSpec.U_sys_nom / config_bat.ModInfo.U_max_mod);  a
    %}

    % Determine number of modules connected in parallel inside the battery system
    
    p_sys_raw = max(ceil(config_bat.SysSpec.I_sys_max / config_bat.ModInfo.I_max_mod), ceil(config_bat.SysSpec.C_sys_min / config_bat.ModInfo.C_mod));
    
    
    % Spatial arrangement of the parallel connection
    
    % We not only can arrange the modules in a simple one dimensional 's*p'
    % grid, but can arrange the modules in all three dimensions. Therefore in
    % this step we create an 'e*pe' connection from the p_sys modules in
    % parallel so that p = e*pe. Those are arranged next to each other so
    % the total number of modules in the battery system is n = e*pe*s.
    
    epe_sys = epe_distribution(p_sys_raw, config_bat.SysSpec.num_higher_p_sys, config_bat.SysSpec.sys_min_e, config_bat.SysSpec.sys_max_e);   
    
    % More parallel connections may be considered for a given configuration
    % --> Interrate through all p's
    
    for jj = 1:1:size(epe_sys,2)   
        
        p_sys = epe_sys(jj).p;  
        
        % For every p several combinations of e-pe may be found. --> We
        % iterate through them as well
        
        for kk = 1:1:size(epe_sys(jj).pe, 2)
            
            config_bat.SysInfo.num_mods_sys = s_sys * p_sys;
            config_bat.SysInfo.num_serial_mods_sys = s_sys;
            config_bat.SysInfo.num_parallel_mods_sys = p_sys;
            config_bat.SysInfo.num_layers_sys = epe_sys(jj).e(kk);
            config_bat.SysInfo.num_parallel_mods_per_layer_sys = epe_sys(jj).pe(kk);
            
            config_bat.SysPara.p_sys = p_sys * config_bat.SysPara.p_mod;
            config_bat.SysPara.s_sys = s_sys * config_bat.SysPara.s_mod;
            config_bat.SysPara.pe_sys = epe_sys(jj).pe(kk) * config_bat.SysPara.pe_mod;
            config_bat.SysPara.e_sys = epe_sys(jj).e(kk) * config_bat.SysPara.e_mod;
            
            config_bat.sys_ID = get_config_ID;  % Assign a unique ID to this config
            
            configs_3_sys_all = append_configs(configs_3_sys_all, config_bat); % Save everything in a cell-array 
        end
    end
end 

% Check if feasible configurations have been found 

check_feasible_configs(configs_3_sys_all,0)


% Save resuls of this step

save(pathes.bat.save_3, 'configs_3_sys_all')  % Save all possible configurations of this run in a MAT-File
%save('BTMS_simulation_results\configs_3_sys_all', 'configs_3_sys_all')  % Save all possible configurations of this run in a MAT-File

clearvars -except config* pathes* % Clear everything instead the array with the configs.


%% Step 4: Basic testing on pack-level

% Similar to step 2 we test the battery packs we created in the last step
% for basic compliance with mass, dimension and energy critera to sort out
% all unfeasible configurations.

clear append_configs         % Clear persistent variable in function
load('configs_3_sys_all');   % Load the feasible configurations from last step if not already in workspace

configs_4_sys_passed = preallocate_configs_4_sys; % Preallocating the cell-array with all configuations that passed the module tests
configs_4_sys_failed = preallocate_configs_4_sys; % Preallocating the cell-array with all configuations that failed the module tests


% Iterate through the configurations

for ii = 1:1:size(configs_3_sys_all, 2)
    
    config_bat = configs_3_sys_all(ii);

    % Exclude Configs that violate the max. dimensions (without consideration of BTMS)

    [config_bat, passed_dim] = sys_test_dimensions_no_BTMS(config_bat);


    % Exclude Configs that violate the energy content criteria

    [config_bat, passed_energy] = sys_test_energy(config_bat);

    % Exclude Configs that violate the max voltage criteria

    [config_bat, passed_voltage] = sys_test_voltage(config_bat);

    % Exclude Configs that can't reach the desired max current criteria

    [config_bat, passed_current] = sys_test_current(config_bat);

    % Create joint structure of criteria

    passed_sys1 = join_passed_structs(passed_dim, passed_energy);
    passed_sys2 = join_passed_structs(passed_voltage, passed_current);
    passed_sys = join_passed_structs(passed_sys1, passed_sys2);


    % Exclude configs that did not pass the tests from further consideration,
    % pass on working configs to the next steps.

    if check_for_failed_tests(passed_sys)  % Test if any test has failed
        configs_4_sys_failed = append_configs(configs_4_sys_failed, config_bat, passed_sys, 'fail');

    else    % Those configurations have passed
        configs_4_sys_passed = append_configs(configs_4_sys_passed, config_bat, passed_sys, 'pass');
    end
end


% Check if feasible configurations have been found 

check_feasible_configs(configs_4_sys_passed,configs_4_sys_failed);

% Save results of this step

save(pathes.bat.save_4_passed, 'configs_4_sys_passed')  % Save all possible configurations of this run in a MAT-File
save(pathes.bat.save_4_failed, 'configs_4_sys_failed')  % Save all failed configurations of this run in a MAT-File
%save('BTMS_simulation_results\configs_4_sys_passed', 'configs_4_sys_passed')  % Save all possible configurations of this run in a MAT-File
%save('BTMS_simulation_results\configs_4_sys_failed', 'configs_4_sys_failed')  % Save all failed configurations of this run in a MAT-File

clearvars -except config* pathes* % Clear everything instead the array with the configs.

try
%% Step 5: Add the BTMS to the modules and system

% Now we have found module and pack configurations that may work from a
% data-sheet based point of view. In this step we take those systems and
% include a BTMS configuration.

clear append_configs            % Clear persistent variable in function
load('configs_4_sys_passed');   % Load the feasible configurations from last step if not already in workspace

configs_5_BTMS_all = preallocate_configs_5_BTMS; % Preallocating the cell-array with all configuations that passed the module tests


% Iterate through the configurations

for ii = 1:1:size(configs_4_sys_passed, 2)
    
    config_bat = configs_4_sys_passed(ii);
    
    % Determine BTMS architecture

    config_bat.BTMSConfig = main_create_BTMS_config(config_bat.BTMSPara, config_bat.SysPara, config_bat.SysInfo, config_bat.ModInfo);
    config_bat.BTMS_ID = config_bat.BTMSPara.name;
    
    % Calculate thermal system properties depending on module and BTMS config
    
%config.SysPara.thermal = calc_thermal_system_properties(config.BTMSConfig, config.SysPara.thermal);
    
    configs_5_BTMS_all = append_configs(configs_5_BTMS_all, config_bat);

end


% Check if feasible configurations have been found 

check_feasible_configs(configs_5_BTMS_all,0);


% Save results of this step

save(pathes.bat.save_5, 'configs_5_BTMS_all')  % Save all possible configurations of this run in a MAT-File
%save('BTMS_simulation_results\configs_5_BTMS_all', 'configs_5_BTMS_all')  % Save all possible configurations of this run in a MAT-File

clearvars -except config* pathes* % Clear everything instead the array with the configs.


%% Step 6: Basic tests of the full system

% Similar to step 4 we test the battery packs for basic compliance with 
% mass, dimension and energy critera - this time with the BTMS included - to
% sort out all unfeasible configurations.

clear append_configs          % Clear persistent variable in function
load('configs_5_BTMS_all');   % Load the feasible configurations from last step if not already in workspace

configs_6_BTMS_passed = preallocate_configs_6_BTMS; % Preallocating the cell-array with all configuations that passed the module tests
configs_6_BTMS_failed = preallocate_configs_6_BTMS; % Preallocating the cell-array with all configuations that failed the module tests


% Iterate through the configurations

for ii = 1:1:size(configs_5_BTMS_all, 2)
    
    config_bat = configs_5_BTMS_all(ii);

    % Exclude Configs that violate the max. dimensions (with consideration of BTMS)

    [config_bat, passed_BTMS] = test_dimensions_BTMS(config_bat);

    % Exclude configs that did not pass the tests from further consideration,
    % pass on working configs to the next steps.

    if check_for_failed_tests(passed_BTMS)  % Test if any test has failed
        configs_6_BTMS_failed = append_configs(configs_6_BTMS_failed, config_bat, passed_BTMS, 'fail');

    else    % Those configurations have passed
        configs_6_BTMS_passed = append_configs(configs_6_BTMS_passed, config_bat, passed_BTMS, 'pass');
    end
end


% Check if feasible configurations have been found 

check_feasible_configs(configs_6_BTMS_passed,configs_6_BTMS_failed);

catch
    err_msg = ['\n########################################################', ...
            '\nNo feasible BTMS solution!\n', '\n########################################################'];
    err_msg = sprintf(err_msg);
    error(err_msg);
end

% Save results of this step

save(pathes.bat.save_6_passed, 'configs_6_BTMS_passed')  % Save all possible configurations of this run in a MAT-File
save(pathes.bat.save_6_failed', 'configs_6_BTMS_failed')  % Save all failed configurations of this run in a MAT-File
%save('BTMS_simulation_results\configs_6_BTMS_passed', 'configs_6_BTMS_passed')  % Save all possible configurations of this run in a MAT-File
%save('BTMS_simulation_results\configs_6_BTMS_failed', 'configs_6_BTMS_failed')  % Save all failed configurations of this run in a MAT-File

clearvars -except config* pathes* % Clear everything instead the array with the configs.

end

%% Next steps

% Refer to 'main_sim_BTMS_2_system_simulation'
