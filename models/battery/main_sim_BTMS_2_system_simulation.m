function [] = main_sim_BTMS_2_system_simulation(battery_GD, config_GD)

%% Info

% Part 2 of the BTMS configuration: Analyze and simulate the systems

% To create the configurations refer to 'main_sim_BTMS_1_system_setup'.

%% GD

% sim_BTMS_2_system_simulation v1.1
% Adapted for use in Global Drive 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                           Versions                                    %
% v1.1 - changed into function                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        known Problems                                 %
% - suboptimal loading of parameters                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization


%%%%%%%%%%%%%%%%%% work in progress
%clearvars
%close all
%clc


%% User Input: Select simulation modes

modes.plot_sys = false;      % Plot the system configuration
modes.sim_module = false;    % Electrical simulation of modules
modes.sim_sys = true;        % Electrical and thermal simulation of the whole system

%Adapted - Global Drive
config_file_mod_sim = 'configs_6_BTMS_passed';   % Select the mat-File with the configs for module simulation
config_file_sys_sim = 'configs_6_BTMS_passed';   % Select the mat-File with the configs for module simulation

% Note: The steps above will be performed for all configs in the provided
% mat-file and therefore may take a long time!

% For further system specifications and simulation parameters go to:
%  - functions_BTMS_sim/sim_module
%  - functions_BTMS_sim/sim_module


%% Plotting of systems

if modes.plot_sys == true

    plot_system_architecture(battery_GD);

    %{
    load(config_file_mod_sim)

    for ii = 1:1:size(configs_6_BTMS_passed, 2)

        config = configs_6_BTMS_passed(ii);
        plot_system_architecture(config);

    end
    
    %}

end


%% Electrical simulation and evaluation on module level

clear sim_module % Clear persistent variable

if modes.sim_module == true
    
    load(config_file_mod_sim)

    for ii = 1:1:size(configs_6_BTMS_passed, 2)

        config = configs_6_BTMS_passed(ii);
        
        results = sim_module(config);
        
        if isempty(results) == false
        
            [results, passed_mod] = check_results_mod(results, config);  

            if check_for_failed_tests(passed_mod)
                
                save_sim_result(results, passed_mod, 'results_7', 'mod', 'failed')
                fprintf('Warning: The module with mod_ID %i failed the electrical tests and should be excluded from further consideration.\n\n', results.mod_ID)
                
            else
                
                save_sim_result(results, passed_mod, 'results_7', 'mod', 'passed')
                
            end
        end
    end
    
    disp('Module simulation finished.')
    clearvars results config t_* SimPara ii passed_mod
    
end


%% Electrical and thermal simulation and evaluation on system level

clear sim_system % Clear persistent variable

if modes.sim_sys == true
    
    load(config_file_sys_sim)

    for ii = 1:1:size(configs_6_BTMS_passed, 2)
        
        results = GD_sim_system(battery_GD, config_GD);

%{        Cleanup
       if isempty(results) == false
        
            [results, passed_sys] = check_results_sys(results, config);  

            if check_for_failed_tests(passed_sys)
                
                save_sim_result(results, passed_sys, 'results_8', 'sys', 'failed')
                fprintf('Warning: The system with sys_ID %i failed the electrical tests and should be excluded from further consideration.\n\n', results.sys_ID)
                
            else
                
                save_sim_result(results, passed_sys, 'results_8', 'sys', 'passed')
                
            end
       end
%}
    end
    
%    disp('System simulation finished.')
%    clearvars results config t_* SimPara ii passed_sys
    
end

%clearvars config_* modes

end