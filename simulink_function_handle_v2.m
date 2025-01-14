function simOut = simulink_function_handle_v2(run_driving_cycle)
%% inputs:
% run_driving_cycle = -1  initialization only (no simulation)
% run_driving_cycle =  0   acceleration cycle only (short simulation)
% run_driving_cycle =  1   driving cycle only (cycle simulation)

% output:
simOut = 0;

%% Close Simulink Model
model = 'BEV_Laengsdynamikmodell_v2';
assignin("base","model",model);
close_system(model, 0);

%% Read vehicle config
config_vehicle = evalin('base','config_vehicle');


%% Battery design variables

    bat_size_factor = 0;
    config_battery = evalin("base", "config_battery");
    bat_celltype = config_battery.Cell_type;
    bat_cellchemistry = config_battery.Cell_chemistry;
    bat_cell_capacity = config_battery.Cell_capacity;
    bat_U_mod_factor = 12;

%% Timing

disp(['Simulation running... ']);
current_time = datetime('now');
time_str = datestr(current_time, 'HH:MM:SS');
disp(time_str);
% t1 = tic; %Timing Initialization


%% Initialize Gearbox

% t_GB = tic; %Timing Gearbox
GB_Mode = evalin("base","GB_Mode"); %Load Gearbox Model type

try
    text = sprintf('Gearbox model: %s',GB_Mode);
    disp(text) %Display model type
    config_gearbox = evalin("base", "config_gearbox"); %Load Gearbox Options
    
    if GB_Mode == "Analytic" %Analytic Gearbox Design
        GB_Topology = evalin("base","GB_Topology");
        if GB_Topology == "achsparallel"
            [GB] = Getriebedesign(config_gearbox.M_max,config_gearbox.ratio,2); %Design "achsparallel"
        elseif GB_Topology == "koaxial"
            [GB] = Getriebedesign_koaxial(config_gearbox.M_max,config_gearbox.ratio); %Design "koaxial"
        elseif GB_Topology == "einstufig"
            [GB] = Getriebedesign(config_gearbox.M_max,config_gearbox.ratio,1); %Design "einstufig"
        end

        assignin('base','GB',GB); %Save Design to Workspace
        GB_Stufen = GB.Stufen; %Extract number of stages
        assignin('base','GB_Stufen',GB_Stufen); %Save number of stages to Workspace

    elseif GB_Mode == "Lookup" %Lookup Gearbox Model
        disp('Loading gearbox lookup...')
        gearboxmap_name = evalin("base","gearboxmap_name");
        load(gearboxmap_name); %Load Lookup-Table
        assignin('base', 'gearbox_map', gearbox_map); %Save Lookup to Workspace
        config_gearbox.iges = gearbox_map.iges;
        assignin('base','config_gearbox',config_gearbox); %Save Design to Workspace
    end

    gear_check = 1; %Initialization Gearbox sucessful
    % disp("Time required for gearbox design: ")
    % toc(t_GB) %Timing Gearbox Design

catch ME
    gear_check = 0;
    display_error_msg(ME);
end


%% Initialize Battery

% tbat = tic; %Timing Battery Design
config_battery = evalin("base", "config_battery");

%line 87-92 + 115-117 only for development: reuse battery model
% Use only, if you want to safe time by using the already built battery model
% try battery_built = evalin("base", "battery_built"); end %try to load "battery_built" variable
% 
% if exist("battery_built") %check, if battery already exists
%    disp("Battery already built, being reused to save time, MAKE SURE TO REMOVE AT THE END!!")
%     battery_check = 1; %skip rebuilt of battery for time saving
% else

    try

    % choose cell
    % config_battery = cell_loader(config_battery, bat_cell_idx);    %#ok<NASGU,NODEF>

    %Battery Optimization %%%
    config_battery = cell_loader(config_battery, bat_cell_capacity, bat_celltype ,bat_cellchemistry);
    try
        config_battery.U_mod_nom = bat_U_mod_factor * 4;
    catch
        config_battery.U_mod_nom = 48;
        warning("Battery: Set U_mod_nom at default value of 48 [V]");
    end

    assignin('base', "config_battery", config_battery)
    % build battery
    battery_check = init_simulation.init_battery(bat_size_factor);
    
    if config_battery.simulation_charge_dynamics == "rc2"
        ModuleType1 = evalin("base","ModuleType1");
        ModuleType1.R3_matCell = ModuleType1.R2_matCell;
        ModuleType1.tau3_matCell = ModuleType1.tau2_matCell;
        ModuleType1.simulation_charge_dynamics = 2;
        assignin('base',"ModuleType1",ModuleType1);
    end

    catch ME
        battery_check = 0;
        display_error_msg(ME);
    end
%     battery_built = true; %See line 83 for details
%     assignin('base',"battery_built",battery_built);
% end

% disp("Time required for battery: ")
% toc(tbat) %Timing Battery Design


%% Initialize Motor

Motor_Mode = evalin("base","Motor_Mode"); %Load Motor Model type

try
    % tmot = tic; %Timing Motor Initialization
    text = sprintf('Motor model: %s',Motor_Mode);
    disp(text) %Display Motor Model type

    %Coefficient for mech. Loss Model: P_FW = k_FW * n^2.5
    %Standard: 410W Mechanical Loss @ 14000rpm
    %Source for Standard Values:
    % [119] Ricciardi, "Development of Virtual Methodology to Evaluate Electric Motor Losses", Master Thesis, Turin, 2023       
    k_FW = 410/(14e3)^2.5;
    assignin('base',"k_FW",k_FW);

    if Motor_Mode == "Lookup" %Lookup Motor Model
        motormap_name = evalin("base","motormap_name");
        load(motormap_name); %Load Lookup-Table for Motor
        MCad_Kennfeld_Vorbereitung; %Prepare Lookup-Table for Simulation
        assignin('base',"motormap",motormap); %Safe Lookup-Table to base Workspace
    end

    if Motor_Mode == "Analytic" %Analytic Motor Model: Design + Analysis
        config_motor = evalin("base","config_motor"); %Load Motor Options
            
            %Prepare Options for handover to MEAPA
            % Motortype 'PMSM' or 'ASM'
            rated.Maschinentyp = config_motor.type;
            % Nominal power P_N [W]
            rated.P_N   = config_motor.power;
            % Nominal speed n_N [U/min]
            rated.n_N   = config_motor.n_nenn;
            % Nominal Voltage U_N [V]
            rated.U_N   = (config_motor.U_nenn-20)/sqrt(2); % = U_LL,eff
            % Number of pole pairs p [-]
            rated.p     = config_motor.polePairs;
            % Rated frequency f_N [-]
            rated.f_N   = (rated.p * rated.n_N) / 60;
            % Number of phases m [-]
            rated.m     = config_motor.phase;
        
        [Entwurf, Analyse] = MEAPA_Skript(rated); %Start MEAPA Tool
        
        %Save Motor Design & Analysis Data to base Workspace
        assignin('base', 'Entwurf', Entwurf);
        assignin('base', 'Analyse', Analyse);
    end
    
    % disp("Time required for motor design & analysis: ")
    % toc(tmot) %Timing Motor Initialization
    motor_check = 1;
    
catch ME
    motor_check = 0;
    display_error_msg(ME);
end


%% Check sucessful Initialization
if gear_check == 0
    disp('Error encountered at gearbox init!')
    disp('Try another combination');
end
if motor_check == 0
    disp('Error encountered at motor init!')
    disp('Try another combination');
end
if battery_check == 0
    disp('Error encountered at battery init!')
    disp('Try another combination');
end

% disp("Time required for initialization: ")
% toc(t1) %Timing Initialization


%% Simulation:

if motor_check && battery_check && gear_check && run_driving_cycle >= 0
    
    if true % Display battery information
        battery = evalin('base', 'battery');
        try
            disp(['Battery size: z=', num2str(battery.SysInfo.dim_z_sys_BTMS*1e3, 4), ' mm'])
        catch 
            disp(['Battery size: z=', num2str(battery.SysInfo.dim_z_sys*1e3, 4), ' mm'])
        end
        disp(['Battery module structure: ', num2str(battery.ModInfo.num_serial_cells_mod),'s', num2str(battery.ModInfo.num_parallel_cells_mod),'p'])
        disp(['Battery system structure: ', num2str(battery.SysInfo.num_serial_mods_sys), 's', num2str(battery.SysInfo.num_parallel_mods_sys),'p'])
        disp(['Battery nominal voltage: U=', num2str(battery.SysInfo.U_nom_sys, 4), ' V'])
        disp(['Battery capacity: E=', num2str(battery.SysInfo.E_sys, 4), ' kWh'])
        disp(['Battery mass: m=', num2str(battery.SysInfo.mass_sys, 4), ' kg'])
    end

    %% Load simulation model
    load_system(model);

    if run_driving_cycle == 0 
        dc_cycle = {'custom'}; %set acceleration cycle
    else
        dc_cycle = {'driving'}; %set driving cycle
    end

    if strcmp(dc_cycle{1}, 'custom') %Fixed battery voltage for acceleration
        acc_cycle_flag = 1;
        disp('Temporarly use fixed battery voltage for acceleration cycle.')
    else
        acc_cycle_flag = 0; %Variable battery voltage for inverter
    end
    assignin('base', 'acc_cycle_flag', acc_cycle_flag); %save to workspace
    

    if strcmp(dc_cycle{1}, 'driving') %Battery overcurrent option
        config_battery.overcurrentfactor_racemode = 1; %for driving cycle
    else
        config_battery.overcurrentfactor_racemode = 2; %for accelerat. cycle
    end
    disp(['Running with battery overcurrent: ', num2str(config_battery.overcurrentfactor_racemode)])


    %% Run simulation
    tsim = tic;
    import_dc_cycle(dc_cycle{1}); %Load relevant driving/acc. cycle
    simOut = sim(model); %Start simulation in Simulink
    disp('Time for simulation:')
    toc(tsim) %Timing Simulink simulation only
    

    %% Gather output
    if strcmp(simOut.SimulationMetadata.ExecutionInfo.StopEvent, 'ModelStop') %if simulation fails
        assignin('base', 'simOut', simOut);
        disp(['Simulation boundaries in ', dc_cycle{1}, ' cycle not fulfilled. Simulation failed.'])

    elseif strcmp(dc_cycle{1}, 'driving') || run_driving_cycle == 0 %if simulation successful
        dc = evalin('base', 'dc'); %Load driving cycle from base workspace
        battery = evalin('base', 'battery'); %Load battery variable from base workspace
        
        sim_out_cell = struct2cell(simOut.sim_out); 
        for i_signals=1:length(sim_out_cell) %Check Simulink signals for NaN values
            if isnan(sum(sim_out_cell{i_signals}.Data))
                disp(['Warning: Found NaN values in column: ', sim_out_cell{i_signals}.Name]);
            end
        end


        %% energy calculations on module and system level

        % calculate driven distance during cycle
        deltaT = 1;
        distance = cumtrapz(dc.speed) * deltaT / 1e3;
        distance = distance(end);

        % SoC difference cycle start to end
        SoC_diff = 100*(simOut.sim_out.bat_SOC.Data(1) - simOut.sim_out.bat_SOC.Data(end));
        disp(['Difference in SoC between start and end: ', num2str(SoC_diff, 3), ' %'])

        % calculate total energy consumption of system (discharge of battery + battery loss)
        total_energy_consumption = -sum(simOut.sim_out.battery_power_in.Data, 'omitnan') * deltaT / 1e3 / 3600 ...
            + sum(simOut.sim_out.battery_powerloss.Data, 'omitnan') * deltaT / 1e3 / 3600;
        energy = total_energy_consumption / distance * 100; %kWh/100km
        disp(['System total energy consumption: ', num2str(total_energy_consumption, 4), ' kWh'])
        disp(['Energy consumption: ', num2str(energy, 3), ' kWh/100km'])
        
        % calculate energy consumption of battery
        total_energy_bat = sum(simOut.sim_out.battery_powerloss.Data, 'omitnan') * deltaT / 1e3 / 3600;
        
        % calculate energy consumption of inverter
        total_energy_inverter_out = sum(simOut.sim_out.electric_power_out.Data, 'omitnan') * deltaT / 1e3 / 3600;
        total_energy_motor_out = sum(simOut.sim_out.terminal_power.Data, 'omitnan') * deltaT / 1e3 / 3600;
        total_energy_inverter = total_energy_inverter_out - total_energy_motor_out;

        % calculate energy consumption of motor
        total_energy_gearbox_out = sum(simOut.sim_out.transmission_power_out.Data, 'omitnan') * deltaT / 1e3 / 3600;
        total_energy_motor = total_energy_motor_out - total_energy_gearbox_out;
        
        % calculate energy consumption of gearbox and vehicle
        total_energy_vehicle = sum(simOut.sim_out.vehicle_power_out.Data, 'omitnan') * deltaT / 1e3 / 3600;
        total_energy_gearbox = total_energy_gearbox_out - total_energy_vehicle;

        % calculate energy consumption of auxiliaries
        drain_aux = config_vehicle.drain_auxiliaries;
        total_energy_aux = drain_aux * (dc.time(end)) / 1e3 / 3600;

        % total_energy_bat = total_energy_bat - total_energy_aux - total_energy_inverter;
        
        disp(['Vehicle total energy consumption: ', num2str(total_energy_vehicle, 4), ' kWh'])
        disp(['Auxiliaries total energy consumption: ', num2str(total_energy_aux, 4), ' kWh'])
        disp(['Gearbox total energy consumption: ', num2str(total_energy_gearbox, 4), ' kWh'])
        disp(['Motor total energy consumption: ', num2str(total_energy_motor, 4), ' kWh'])
        disp(['Inverter total energy consumption: ', num2str(total_energy_inverter, 4), ' kWh'])
        disp(['Battery total energy consumption: ', num2str(total_energy_bat, 4), ' kWh'])
        
        % calculate recuperation at battery level
        recup = simOut.sim_out.battery_power_out.Data .* deltaT;
        recup = recup(recup > 0); %only recuperation power >0
        recup = sum(recup) / 1e3 / 3600;
        disp(['Recuperation energy: ', num2str(recup, 3), ' kWh'])

        % calculate average power loss in kW
        power_bat = total_energy_bat / (dc.time(end)) * 3600;
        power_inverter = total_energy_inverter / (dc.time(end)) * 3600;
        power_mot = total_energy_motor / (dc.time(end)) * 3600;
        power_gear = total_energy_gearbox / (dc.time(end)) * 3600;
        power_vehicle = total_energy_vehicle / (dc.time(end)) * 3600;

        disp(['Vehicle average power loss: ', num2str(power_vehicle, 3), ' kW'])
        disp(['Auxiliaries average power loss: ', num2str(drain_aux / 1e3, 3), ' kW'])
        disp(['Gearbox average power loss: ', num2str(power_gear, 3), ' kW'])
        disp(['Motor average power loss: ', num2str(power_mot, 4), ' kW'])
        disp(['Inverter average power loss: ', num2str(power_inverter, 3), ' kW'])
        disp(['Battery average power loss: ', num2str(power_bat, 3), ' kW'])


        % save simulation output to workspace
        assignin('base', 'simOut', simOut)
    end
end

end


function display_error_msg(ME)
disp('Error encountered!')
disp(ME.message)
for e_i = 1:length(ME.stack)
    disp(ME.stack(e_i))
end
end
