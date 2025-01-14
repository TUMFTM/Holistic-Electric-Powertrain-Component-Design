function [] = battery_builder(battery_GD, config, pathes, visualization)
    arguments
        battery_GD struct = struct();
        config struct = struct();
        pathes struct = struct();
        visualization struct = table2struct(table(0,0,0,0,0, ...
        'VariableNames',["Cell" "Module_P" "Module_S" "Module_A" "Pack"]));
    end


%% battery_builder v0.12
% This script builds the battery according to the given Specifications
% part of Global Drive 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                           Versions                                    %
% v0.1 - loading cell data - prismatic only                             %
% v0.2 - added module - ParallelAssembly                                %
% v0.3 - added Module - SerialAssembly                                  %
% v0.4 - added Script Checker for better testing                        %
% v0.5 - added Module Assembly                                          %
% v0.6 - added Pack                                                     %
% v0.7 - Generate Simscape Battery library model                        %
% v0.8 - added passiv cell balancing                                    %
% v0.9 - added Thermal Effects - beta                                   %
% v0.10 - added charge dynamics                                         %
% v0.11 - preparations for merge                                        %
% v0.12 - added compability with Pre-Parameterized Cells                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        known Problems                                 %
% - Module Assembly and Pack weight and volume are off                  %
% - geometry of module arrangement is off  -> implemtent z dimension    %
% - Pack builder not finished jet - focused on testing first            %
% - cleanup needed                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scriptchecker = 0;

%% init
import simscape.battery.builder.*

%% geometry
batteryCell = Cell();


if ismember(battery_GD.BatPara.cell_type, ["Cyl", "Cylindrical"])
    cellGeometry = CylindricalGeometry();
    cellGeometry.Radius = simscape.Value(battery_GD.BatPara.physical.dim_x,"m");
    cellGeometry.Height = simscape.Value(battery_GD.BatPara.physical.dim_z, "m");

    %cellGeometry = CylindricalGeometry();
    %cellGeometry.Radius = simscape.Value(0.0105, "m");
    %cellGeometry.Height = simscape.Value(0.07, "m");

    %{
    batteryCell.CellModelOptions.BlockParameters.thermal_port = "model";

    batteryParallelAssembly = ParallelAssembly(Cell = batteryCell,...
    NumParallelCells = 4, ...
    Rows = 4, ...
    Topology = "Square", ...
    ModelResolution = "Detailed");

    batteryModule = Module(ParallelAssembly = batteryParallelAssembly,...
    NumSeriesAssemblies = 13, ...
    InterParallelAssemblyGap = simscape.Value(0.005,"m"), ...
    ModelResolution = "Detailed", ...
    AmbientThermalPath = "CellBasedThermalResistance");

    disp(batteryModule.NumModels);

    disp(batteryModule.NumModels);

    nexttile(tl)
    batteryModuleChart1 = BatteryChart(Parent = tl, Battery = batteryModule);
    nexttile(tl)
    batteryModuleChart2 = BatteryChart(Parent = tl, Battery = batteryModule, SimulationStrategyVisible = "On");
    %}

elseif battery_GD.BatPara.cell_type == "Pouch"
    cellGeometry = PouchGeometry();
    cellGeometry.Thickness = simscape.Value(battery_GD.BatPara.physical.dim_x,"m");
    cellGeometry.Length = simscape.Value(battery_GD.BatPara.physical.dim_y, "m");
    cellGeometry.Height = simscape.Value(battery_GD.BatPara.physical.dim_z, "m");
    cellGeometry.TabLocation = battery_GD.BatPara.physical.TabLocation;
    cellGeometry.TabWidth = simscape.Value(battery_GD.BatPara.physical.TabWidth, "m");
    cellGeometry.TabHeight = simscape.Value(battery_GD.BatPara.physical.TabHeight, "m");

elseif battery_GD.BatPara.cell_type == "Prismatic"
    cellGeometry = PrismaticGeometry();
    cellGeometry.Thickness = simscape.Value(battery_GD.BatPara.physical.dim_x,"m");
    cellGeometry.Length = simscape.Value(battery_GD.BatPara.physical.dim_y, "m");
    cellGeometry.Height = simscape.Value(battery_GD.BatPara.physical.dim_z, "m");

else
    error("unsupported celltype/n check configuration file")
    
end

batteryCell.Geometry = cellGeometry;

batteryCell.Mass = simscape.Value(battery_GD.BatPara.physical.m,"kg");
batteryCell.CellModelOptions.BlockParameters.thermal_port = "model";    % Thermal modeling
batteryCell.CellModelOptions.BlockParameters.T_dependence = "yes";      % Thermal modeling

if isempty(strfind(config.simulation_charge_dynamics,"rc"))
    error("unsupported charge dynamics - pls check configuration file")
end

if config.simulation_charge_dynamics == "rc0"
    config.simulation_charge_dynamics = "off";   
end
batteryCell.CellModelOptions.BlockParameters.prm_dyn = config.simulation_charge_dynamics; % Simulate Charge Dynamics [off, rc1, rc2, rc3, rc4, rc5] 

% if config.manufacturer_parameterization == 0
%     % batteryCell.ParameterizationManufacturer = config.manufacturer;
%     % batteryCell.ParameterizationPartNumber = config.cellname;
%     % batteryCell.applyCellDataFromPart;
% 
%     % if battery_GD.BatPara.cell_type == "Cyl"
%     %     batteryCell.Geometry = CylindricalGeometry();
%     % elseif battery_GD.BatPara.cell_type == "Pouch"
%     %     batteryCell.Geometry = PouchGeometry();  
%     % elseif battery_GD.BatPara.cell_type == "Prismatic"
%     %     batteryCell.Geometry = PrismaticGeometry();
%     % end
% 
% end


if visualization.Cell == 1
    f = uifigure("Color", "white");
    cellChart = BatteryChart(Parent = f, Battery = batteryCell);
    title(cellChart, "Pouch Cell")
end

% Create Thermic 1-D Model of battery - not used now
%batteryCell.CellModelOptions.BlockParameters.thermal_port = "model";

%% module - ParallelAssembly
parallelAssembly = ParallelAssembly(...
    ModelResolution=config.simulation_cell_grouping,...
    NumParallelCells = battery_GD.ModInfo.num_parallel_cells_mod, ...
    Cell = batteryCell, ...
    InterCellGap = simscape.Value(battery_GD.BatPara.physical.dim_x * battery_GD.SysSpec.sf_dim_mod * ...
    (battery_GD.ModInfo.num_parallel_cells_mod + 1) / battery_GD.ModInfo.num_parallel_cells_mod, "m"));

if ismember(battery_GD.BatPara.cell_type, ["Cyl", "Cylindrical"])
    parallelAssembly.Topology = "Square";
else
    parallelAssembly.Topology = "SingleStack";
end

disp("battery_builder: Successfully build parallelAssembly")

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Notes: %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
% ModelResolution
% "Lumped" = 1 elcetrical modell and only thermal mass or 
% "Detailed" = p electrical and thermal models
%%%
% NumParallelCells:
% how many cells are electrically connected in a parallel configuration 
% -> no intel about the topology of the assembly
%%%
% Rows =
% Number of rows of the parallel assembly relative to the stacking axis
%%%
% Topolgy =
% Geometrical arrangement of the cells relative to the cell format,
%%%
% MassFactor
% Factor for aditional non-cell-related mass (tabs, collector, ...)
%%%
% InterParallelAssemblyGap = Gap between parallel assembly parts
% Using the SF (=Safetyfactor) corrected for the "in between cells" error
%%%%%%%%%%%%%%%%%%%%%%

if visualization.Module_P == 1
    f = uifigure("Color", "white");
    parallelAssemblyChart = BatteryChart(Parent = f, Battery = parallelAssembly);
    title(parallelAssemblyChart, "Parallel Assembly Chart")
end

%% Module - SerialAssembly
module = Module(...
    ParallelAssembly = parallelAssembly, ...
    MassFactor = 1 + battery_GD.SysSpec.sf_mass_mod,... %adding mass for passive parts
    ModelResolution = config.simulation_cell_grouping,...
    NumSeriesAssemblies = battery_GD.ModInfo.num_serial_cells_mod, ...
    InterParallelAssemblyGap = simscape.Value(battery_GD.BatPara.physical.dim_x * battery_GD.SysSpec.sf_dim_mod *...
    (battery_GD.ModInfo.num_serial_cells_mod + 1) / battery_GD.ModInfo.num_serial_cells_mod, "m"));

disp("battery_builder: Successfully build module")

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Notes: %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%%%
% NumSeriesAssemblies:
% how many cells are electrically connected in a serial configuration 
% -> no intel about the topology of the assembly
%%%
% MassFactor
% Factor for aditional non-cell-related mass (tabs, collector, ...)
%%%
% InterParallelAssemblyGap = Gap between parallel assembly parts
% Using the SF (=Safetyfactor) corrected for the "in between cells" error
%%%%%%%%%%%%%%%%%%%%%%

if visualization.Module_S == 1
    f = uifigure("Color", "white");
    moduleChart = BatteryChart(Parent = f, Battery = module);
    title(moduleChart, "Module Chart")
end

%% Module Assembly
% Assembly of in series connected modules
% "By default, the ModuleAssembly object electrically connects the modules in series."

moduleAssembly = ModuleAssembly(...
    Module = repmat(module,1,battery_GD.SysInfo.num_serial_mods_sys), ...
    CircuitConnection = "Series",...
    NumLevel = battery_GD.SysInfo.num_layers_sys,...
    InterModuleGap = simscape.Value(0.1, "m")); % needs addaption
    % Module = repmat(module,1,battery_R.SysInfo.num_serial_mods_sys), ...

disp("battery_builder: Successfully build moduleAssembly")

if visualization.Module_A == 1 
    f = uifigure("Color", "white");
    moduleAssemblyChart = BatteryChart(Parent = f, Battery = moduleAssembly);
    title(moduleAssemblyChart, "Module Assembly Chart")
end

%% Pack

batteryPack = Pack(...
    ModuleAssembly = [repmat(moduleAssembly,1,battery_GD.SysInfo.num_parallel_mods_sys)], ...
    CircuitConnection = "Parallel",...
    StackingAxis = "X",...
    MassFactor = 1 + battery_GD.SysSpec.sf_mass_mod,... %adding mass for passive parts
    InterModuleAssemblyGap = simscape.Value(0.005, "m"));

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Notes: %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%%%
% CircuitConnection = "Series" (default) 
%%%
% MassFactor
% Factor for aditional non-cell-related mass (tabs, collector, ...)
%%%
% InterParallelAssemblyGap = Gap between parallel assembly parts
% Using the SF (=Safetyfactor) corrected for the "in between cells" error
%%%%%%%%%%%%%%%%%%%%%%

% Setting Balancing Strategy
%batteryPack.BalancingStrategy = "";        %= no balancing
batteryPack.BalancingStrategy = "Passive";   %= passive balancing
moduleAssembly.BalancingStrategy = "Passive";   % Retest if this is really needed
parallelAssembly.BalancingStrategy = "Passive"; % Retest if this is really needed

% Thermal
batteryPack.AmbientThermalPath = "CellBasedThermalResistance";
batteryPack.CoolantThermalPath = "CellBasedThermalResistance";

moduleAssembly.AmbientThermalPath = "CellBasedThermalResistance";
moduleAssembly.CoolantThermalPath = "CellBasedThermalResistance";

module.AmbientThermalPath = "CellBasedThermalResistance";
module.CoolantThermalPath = "CellBasedThermalResistance";
module.CoolingPlate = "Bottom";

parallelAssembly.AmbientThermalPath = "CellBasedThermalResistance";
parallelAssembly.CoolantThermalPath = "CellBasedThermalResistance";

disp("battery_builder: Successfully build pack")

% Visualization
if visualization.Pack == 1
    f = uifigure("Color", "white");
    packChart = BatteryChart(Parent = f, Battery = batteryPack);
    title(packChart, "Pack Chart")
end

%% Script Checker
if scriptchecker == 1
    fprintf("\nScript Checker activated\n")
    
    % module mass check
    fprintf("module mass diff = %d\n\n", (battery_GD.ModInfo.mass_mod - value(module.CumulativeMass)));

    % module volume
    mod_Volume_noBTMS = battery_GD.ModInfo.dim_x_mod * battery_GD.ModInfo.dim_y_mod* battery_GD.ModInfo.dim_z_mod;
    
    fprintf("module volume no BTMS diff = %d \n", (mod_Volume_noBTMS - value(module.PackagingVolume)));
    
    % not looking at BTMS volume as in Reiter
    % mod_Volume_BTMS = battery_R.ModInfo.dim_x_mod_BTMS * battery_R.ModInfo.dim_y_mod_BTMS * battery_R.ModInfo.dim_z_mod_BTMS;
    % fprintf("module volume with BTMS diff = %d \n", (mod_Volume_BTMS -
    % value(module.PackagingVolume)));

    % number of cells in module
    if battery_GD.ModInfo.num_cells_mod == parallelAssembly.NumParallelCells * module.NumSeriesAssemblies
        fprintf("\nNumber of module cells unchanged\n")
    else
        fprintf("\n-!-!-!-!-\nNumber of module cells DOES NOT FIT\n-!-!-!-!-\n")
    end

    
    % pack mass check
    fprintf("module mass diff = %d\n\n", (battery_GD.SysInfo.mass_sys - value(batteryPack.CumulativeMass)));

    % pack volume
    % mod_Volume_noBTMS = battery_R.SysInfo.dim_x_sys * battery_R.SysInfo.dim_y_sys * battery_R.ModInfo.dim_z_sys; % not adapted jet

    % number of cells in module
    % if battery_R.SysInfo
    % end    

end

%% Generate Simscape Battery library model

% Create new folder if its missing
if ~isfolder(pathes.bat.Batterymodel)
    mkdir(pathes.bat.Batterymodel)
    addpath(pathes.bat.Batterymodel)
end

% clear folder of previous model
if isfolder(pathes.bat.Batterymodel + "\+" + config.libraryname)
    rmdir(pathes.bat.Batterymodel + "\+" + config.libraryname,'s')
end

if isfolder(pathes.bat.Batterymodel + "\+" + config.libraryname + "LumpingAdapters")
        rmdir(pathes.bat.Batterymodel + "\+" + config.libraryname + "LumpingAdapters",'s')
end

if ~numel(dir(pathes.bat.Batterymodel))<=2
    delete(pathes.bat.Batterymodel + "\*")
end

% generate new battery pack model
%buildBattery(batteryPack,LibraryName = libraryname,Directory = pathes.bat.Batterymodel);

% build with masked parameters
buildBattery(batteryPack,LibraryName = config.libraryname,Directory = pathes.bat.Batterymodel,...
    MaskInitialTargets = "VariableNames", MaskParameters = "VariableNames");

% buildBattery(batteryPack,LibraryName = config.libraryname,Directory = pathes.bat.Batterymodel);

disp("battery_builder: Successfully generated new Simulink library")

end
