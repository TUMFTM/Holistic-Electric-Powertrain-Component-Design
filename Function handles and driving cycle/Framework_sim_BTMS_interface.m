function [battery_GD] = Framework_sim_BTMS_interface()
    %% Framework_sim_BTMS_interface v0.15
    % This interface for sim_BTMS was created as part of Global Drive 2023
    %
    % Purpose:
    %   Interface for battery thermal management simulation (sim_BTMS).
    % 
    % V3: This function has been greatly restructured for the v3. Many of
    % its original contents can now be found in subfunctions run in the
    % powertrain design Framework.

    %% Test Parameters
    disp('The function has been called!');

    % Skipper for faster testing
    testing_skipper = 0;
    config_battery = evalin('base', 'config_battery');
    battery_GD = evalin("base", "battery_GD");
    BatPara = evalin("base", "BatPara");

    % Skip to specific step for testing purposes
    if testing_skipper > 0
        fprintf("!!!\n!!!\n\nSkipping sim_BTMS_1 directly to step %d\n\n!!!\n!!!\n", testing_skipper);
        configs_6_BTMS_passed = "configs_6_BTMS_passed.mat";
        load(configs_6_BTMS_passed);
        battery_GD.config = configs_elector(battery_GD, configs_6_BTMS_passed);
    end

    %% Step 1: Initialization
    import simscape.battery.builder.*;

    visualization.Cell = 0;
    visualization.Module_P = 0;
    visualization.Module_S = 0;
    visualization.Module_A = 0;
    visualization.Pack = 0;

    % System setup for paths
    pathes.bat = jsondecode(fileread(fullfile("input_and_parameters", "00_programm_pathes", "sim_BTMS_pathes.json")));
    pathes.bat = pathbuilder(pathes.bat, fullfile("models", "battery"));

    %% Step 4: Simscape Battery
    if config_battery.simulation_simscape == 1 && testing_skipper < 4
        % Build battery pack and initialize it
        battery_builder(battery_GD, config_battery, pathes, visualization);
        battery_GD = battery_initializer(battery_GD, config_battery, pathes);

        % Replace existing battery pack and add loss probe
        battery_PackReplacement = battery_PackReplacer();
        battery_LossProbe = battery_LossProbeAdder(battery_GD);

        % Placeholder for running battery simulation
    end

    %% Step 5: Reiter Simulation
    if config_battery.simulation_reiter == 1 && testing_skipper < 4
        main_sim_BTMS_2_system_simulation(battery_GD, config_battery);
    end

    %% Helper Functions
    function [BatPara] = get_BatPara(BatPara, config)
        % Load battery parameters based on configuration
        try
            run(config.cellname);
        catch
            % Fallback to loading from predefined path
            load(fullfile("models", "battery", "input_and_parameters", "01_cell_data", config.cellname, config.cellname + ".mat"));
        end
    end

    % Function placeholder (commented for now)
    % function [] = sim_BTMS_1_caller(BatPara, config, input_configs, pathes, SysSpec)
    %     run main_sim_BTMS_1_system_setup
    % end
end