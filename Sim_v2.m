% This is a script to setup and run a simulation in v2

% optional: clearing the workspace
% clear all
% close all

tgesamt = tic; %Timing for Sim_v2 script

%% Initialize configuration
% Add subfolders to project path
addpath(genpath("config"));
addpath(genpath("data"));
addpath(genpath("models"));

% Import simulation model configuration file
load("SharedConfig.mat");


%% ##### !!User changeable Setup here!! ##################################

% Define if driving cycle should be run
% -1: initialize only
%  0: Accelerationcycle only
%  1: Driving Cycle (only)
run_driving_cycle = 1;


%% Setup for Vehicle & Driveing Cycle:

% Select vehicle file:
% vehicle_config_file = 'config/conf_vehicle_new.json';
vehicle_config_file = 'config/conf_vehicle_ID3.json';
% vehicle_config_file = 'config/conf_vehicle_Tesla3.json';

% Import vehicle configuration (don't change):
disp(['Loading vehicle config: ', vehicle_config_file]);
config_vehicle = read_json_config(vehicle_config_file);

% Set specific details for Accelerationcycle here:
% config_vehicle.corner_speed_acc_cycle = 100; %km/h
% config_vehicle.zero_to_corner_speed_in_s = 7.3; %s
% config_vehicle.maximum_speed = 160; %km/h

% Set specific details for Driving Cycle here:
% 'WLTP_class_3', 'USA_FTP72', 'USA_FTP72', 'EUROPE_NDEC', 'JAPAN_10_MODE',...
% For all available cycles check data/DrivingCycles folder
% config_vehicle.driving_cycle = 'WLTP_class_3';

% Change specific vehicle model details here:
% config_vehicle.mass_total = 1976; %kg
% config_vehicle.cw_air = 0.2627;
% config_vehicle.frontal_area = 2.3904; %m^2
% config_vehicle.drain_auxiliaries = 280; %W
% config_vehicle.rot_mass_factor = 1.07;
% config_vehicle.tire_radius = 0.335; %m
% config_vehicle.tire_roll_resistance = 0.009;


%% Setup for Gearbox:
GB_Mode = "Analytic"; %Gearbox Model choice: 'Analytic' 'Lookup' 'Simple' 'Const'

if GB_Mode == "Analytic" %Setup for analytic model:
    GB_Topology = 'achsparallel'; % Choose topology: 'einstufig','achsparallel' oder 'koaxial'
    config_gearbox.ratio = 11.53; % Gearbox total ratio
    config_gearbox.M_max = 310; %Nm; Maximum input torque
    config_gearbox.power_oilpump = 500; %Watt; maximum power oilpump

elseif GB_Mode == "Lookup" %Setup for lookup model:
    gearboxmap_name = 'transmission_efficiency_map_TP3_11_53.mat'; %Input lookup name
    config_gearbox.power_oilpump = 500; %Watt; maximum power oilpump

elseif GB_Mode == "Simple" %Setup for simple model:
    config_gearbox.ratio = 11.53; % Gearbox total ratio
    config_gearbox.k_n1 = 1.7e-4; % coefficient for speed-proportional losses (linear)
    config_gearbox.k_n2 = 3.0e-9; % coefficient for speed-proportional losses (quadratic)
    config_gearbox.k_M1 = 1.7e-3; % coefficient for torque-proportional losses (lin.)
    config_gearbox.k_M2 = 1.8e-7; % coefficient for torque-proportional losses (quad.)

elseif GB_Mode == "Const" %Setup for constant model:
    config_gearbox.ratio = 11.53; % Gearbox total ratio
    config_gearbox.efficiency = 1; % efficiency 0 < eta <= 1
end


%% Setup for Motor:
Motor_Mode = "Lookup"; %Motor Model Choice: 'Analytic' 'Lookup' 'Const'

if Motor_Mode == "Analytic" %Setup for analytic model:
    config_motor.power = 70e3; %W; Nominal Power
    config_motor.n_nenn = 4600; %1/min; Nominal Speed
    config_motor.type = 'PMSM'; % Motor type: 'PMSM' oder 'ASM'
    config_motor.U_nenn = config_battery.U_sys_nom; %V; Nominal Voltage
    config_motor.polePairs = 4; % Pole pair number
    config_motor.phase = 3; % Phase number (MEAPA needs 3)
    config_motor.M_nenn = config_motor.power/(2*pi*config_motor.n_nenn/60); %Nm; Nominal Torque

elseif Motor_Mode == "Lookup" %Setup for lookup model:
    motormap_name = 'motor_map_PSM_ID3_Deininger.mat'; %Input lookup name
    config_motor.type = 'PMSM'; % Motor type: 'PMSM' oder 'ASM'
    config_motor.skalieren = false; %Activate Scale function: true/false
    config_motor.k_skal_n = 1.1; % coefficient for speed
    config_motor.k_skal_M = 1.1; % coefficient for torque

elseif Motor_Mode == "Const" %Setup for constant model:
    config_motor.efficiency = 0.94; %efficiency 0 < eta <= 1
    config_motor.type = 'PMSM'; % Motor type: 'PMSM' oder 'ASM'
end


%% Setup for Inverter:
Inverter_Mode = "Analytic"; %"Analytic" or "Simple"

% Load standard configuration
config_Inverter = read_json_config("config/conf_Inverter.json");

% Both options only valid for analytic model:
config_vehicle.inverter_type = 'MOSFET'; % Inverter Type: 'IGBT' or 'MOSFET'
config_inverter.fswitch = 5000; % Inverter Switching Frequncy


%% Setup for Battery:
Battery_Mode = "AnalyticRC3"; 
%"AnalyticRC3" for analytic with RC3; "AnalyticR0" simplified with R0; "Simscape" for Simscape

% Load standard configuration
config_battery = read_json_config("config/conf_battery.json");

config_battery.I_sys_max = 250; %A: Max. current on system level
config_battery.U_sys_nom = 400; %V: Nominal system voltage
config_battery.E_sys_min = 58; %kWh: Stored energy
config_battery.Cell_type = 2; %Cell Type: 1=cyl; 2=pouch; 3=pris
config_battery.Cell_chemistry = 3; %Chemistry: 1=NMC111; 2=NMC622; 3=NMC721; 4=NMC811; 5=LFP
config_battery.Cell_capacity = 30; %Ah: Cell capacity



% ###########################################################################
% ######## Stop, don't change anything below this line! #####################
% ###########################################################################

%% Battery Pre-Initialization

% Additional Battery things:
config_battery.libraryname = "GD_Batterymodel_NSGAII_1";   %Name of created battery library

load("available_cells.mat")
config_battery.num_available_cells = size(available_cells,1);   % Get number of available cells
disp(['Found n battery cells: ', num2str(config_battery.num_available_cells)])
clear available_cells

% Bat for EM
if config_battery.U_sys_nom == 400 %Voltage for first design [V]
    config_battery.U_sys_design = 360; 
elseif config_battery.U_sys_nom == 800
    config_battery.U_sys_design = 648;
else
    disp("Unknown value for U_sys_nom. Use either 400 or 800 (V).")
    keyboard
end


%% Inverter Initialization

Mode_IM_IGBT = Simulink.Variant('Mode == 1');
Mode_IM_MOSFET = Simulink.Variant('Mode == 2');
Mode_PSM_IGBT = Simulink.Variant('Mode == 3');
Mode_PSM_MOSFET = Simulink.Variant('Mode == 4');
if config_motor.type == "ASM" && config_vehicle.inverter_type == "IGBT"
    Mode = 1;
    warning('off', 'Simulink:Commands:ParamUnknown');
elseif config_motor.type == "ASM" && config_vehicle.inverter_type == "MOSFET"
    Mode = 2;
    warning('off', 'Simulink:Commands:ParamUnknown');
elseif config_motor.type == "PMSM" && config_vehicle.inverter_type == "IGBT"
    Mode = 3;
    warning('off', 'Simulink:Commands:ParamUnknown');
elseif config_motor.type == "PMSM" && config_vehicle.inverter_type == "MOSFET"
    Mode = 4;
    warning('off', 'Simulink:Commands:ParamUnknown');
end



%% Run Simulation


%Start of Simulation
out = simulink_function_handle_v2(run_driving_cycle);


disp("Time for Initialization and Simulation: ")
toc(tgesamt) %Timing for Sim_v2_1 script



%% helper functions
function config_x = read_json_config(file_path)
    % Read json formatted table of config data and import it to Matlab workspace
    fid = fopen(file_path, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    config_x = jsondecode(str);
end