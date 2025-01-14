function [battery_GD] = battery_initializer(battery_GD, config, pathes)
    arguments
        battery_GD struct = struct();
        config struct = struct();
        pathes struct = struct();
    end

%% battery_initializer v0.12
% This script realises the initialisation of the battery
% part of Global Drive 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                           Versions                                    %
% v0.1 - loading parameters as choosen inside the battery builder       %
% v0.2 - use pOCV and HPPC data                                         %
% v0.3 - generic path usage                                             %
% v0.4 - added SOV deviation                                            %
% v0.5 - added saving of results                                        %
% v0.6 - fixed parametrisation names                                    %
% v0.7 - added battery_GD.SysInfo.dim_surface for heat convection       %
% v.08 - fixed Pack Capacity + prep variable current for EM             %
% v0.9 - added possibility to use PreParametrized simscape cells        %
% v0.10 - added support for lumped ParallelAssemblys/Modules            %
% v0.11 - aktivated R_loss                                              %
% v0.12 - deaktivated R_loss                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        known Problems                                 %
% - Preliminary Solution for ParallelAssemblyType1                      %
% - SOC deviator not tested jet                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load builder parameterset
%run(config.libraryname+"_param")

%% Initialise
% h = help variable
h.n_RC = config.simulation_charge_dynamics(3); % get RC-model from config  

% get from config oder battery_GD

%% load cellparameters

if config.manufacturer_parameterization == 0
    % OCV
    load(pathes.bat.Celldata+"/"+config.cellname+"/pOCV.mat");
    if pOCV.SOC_vec(end) > 1
    pOCV.SOC_vec = pOCV.SOC_vec'/100;
    end
    
    % xRC
    x_RC = load(pathes.bat.Celldata+"/"+config.cellname+"/"+h.n_RC+"RC/"+h.n_RC+"RC.mat");
    x_RC =  x_RC.("results_" + h.n_RC + "RC");
    
    %R_loss - For Preliminary loss calculation
    % try
    %     R_loss = load(pathes.bat.Celldata+"/"+config.cellname+"/0RC/0RC.mat");
    % 
    %     R_loss =  R_loss.("results_0RC");
    % catch
    %     fprintf("No R0-Parameters for R_loss found")
    %     R_loss = 0;
    % end

    if isfield(x_RC, 'R0_mat')
        for d_ch = ["ch", "dch"]
            x_RC.("R0_mat_"+d_ch) = x_RC.("R0_mat");
            for i = 1:str2double(h.n_RC)
                x_RC.("R"+string(i)+"_mat_"+d_ch) = x_RC.("R"+string(i)+"_mat");
                x_RC.("tau"+string(i)+"_mat_"+d_ch) = x_RC.("tau"+string(i)+"_mat");
            end
        end 
    end

    if size(x_RC.T_vec,2) == 1
        x_RC.T_vec = x_RC.T_vec + [-25 0 25];    
        for d_ch = ["ch", "dch"]
            x_RC.("R0_mat_"+d_ch) = repmat(x_RC.("R0_mat_"+d_ch),1,numel(x_RC.T_vec));
            for i = 1:str2double(h.n_RC)
                if size(x_RC.("R"+string(i)+"_mat_"+d_ch)) ~= size(x_RC.T_vec,2)
                    x_RC.("R"+string(i)+"_mat_"+d_ch) = repmat(x_RC.("R"+string(i)+"_mat_"+d_ch),1,numel(x_RC.T_vec));
                    x_RC.("tau"+string(i)+"_mat_"+d_ch) = repmat(x_RC.("tau"+string(i)+"_mat_"+d_ch),1,numel(x_RC.T_vec));
                end
            end
       end
    end
    
else
    pOCV.SOC_vec = (battery_GD.BatPara.electrical.SOC)';
    pOCV.U_SOC = (battery_GD.BatPara.electrical.OCV)';

    x_RC.T_vec = [0 25 50 ] + 273.15;                   % Temperature in Kelvin
    x_RC.R0_mat_dch = repmat(battery_GD.BatPara.electrical.R0,size(x_RC.T_vec,2),1).';
    x_RC.R0_mat_ch = x_RC.R0_mat_dch;
    
    % R_loss = x_RC.R0_mat_dch;
end

% Parameters of created Battery Library
if h.n_RC > 0
    h.PAT1_fieldnames_RC = [];
    for i = 1:str2double(h.n_RC)
        h.PAT1_fieldnames_RC = [h.PAT1_fieldnames_RC, "R"+i+"_matCell","tau"+i+"_matCell"];
    end
end

h.ParallelAssemblyType1_fieldnames =["SOC_vecCell","T_vecCell","V0_matCell",...
"V_rangeCell","R0_matCell","AHCell",h.PAT1_fieldnames_RC,...
"thermal_massCell","CoolantResistance","AmbientResistance",...
"CellBalancingClosedResistance","CellBalancingOpenConductance",...
"CellBalancingThreshold","CellBalancingResistance"];

for d_ch = ["ch", "dch"] % charge or discharge
    d.(d_ch).direction = d_ch;
    %% ModuleType1
    
    % Values from cell input
    d.(d_ch).ModuleType1.SOC_vecCell = pOCV.SOC_vec';   % Vector of state-of-charge values, SOC
    d.(d_ch).ModuleType1.T_vecCell = x_RC.T_vec;            % Temperatur of HPPC fit
    %ModuleType1.T_vec = getfield(("results_" + h.n_RC + "RC"),T_vec);
    d.(d_ch).ModuleType1.V0_matCell = repmat(pOCV.U_SOC,1,size(d.(d_ch).ModuleType1.T_vecCell,2)); % As we only got an OCV for 1 Temp
    
    if config.manufacturer_parameterization == 0
        if x_RC.SOC_vec(1) > 0
            d.(d_ch).ModuleType1.R0_matCell =[x_RC.("R0_mat_"+d_ch)(1,:); ...
            x_RC.("R0_mat_"+d_ch)];   % Terminal resistance, R0(SOC,T), Ohm
            for i = 1:str2double(h.n_RC)
                d.(d_ch).ModuleType1.("R"+string(i)+"_matCell") = ...
                [x_RC.("R"+string(i)+"_mat_"+d_ch)(1,:); ...
                x_RC.("R"+string(i)+"_mat_"+d_ch)];
                d.(d_ch).ModuleType1.("tau"+string(i)+"_matCell") = ...
                [x_RC.("tau"+string(i)+"_mat_"+d_ch)(1,:); ...
                x_RC.("tau"+string(i)+"_mat_"+d_ch)];
            end
        end

        if x_RC.SOC_vec(end) < 100
            d.(d_ch).ModuleType1.R0_matCell =[d.(d_ch).ModuleType1.R0_matCell; ...
            d.(d_ch).ModuleType1.R0_matCell(end,:)];   % Terminal resistance, R0(SOC,T), Ohm    
            for i = 1:str2double(h.n_RC)
                d.(d_ch).ModuleType1.("R"+string(i)+"_matCell") = ...
                [d.(d_ch).ModuleType1.("R"+string(i)+"_matCell");...
                d.(d_ch).ModuleType1.("R"+string(i)+"_matCell")(end,:)];
                d.(d_ch).ModuleType1.("tau"+string(i)+"_matCell") = ...
                [d.(d_ch).ModuleType1.("tau"+string(i)+"_matCell"); ...
                d.(d_ch).ModuleType1.("tau"+string(i)+"_matCell")(end,:)];
            end
        end
        
        % if isstruct(R_loss)
        %     d.(d_ch).R_loss =[R_loss.("R0_mat_"+d_ch)(1,:);R_loss.("R0_mat_"+d_ch);...
        %     R_loss.("R0_mat_"+d_ch)(end,:)];
        % else
        %     d.(d_ch).R_loss = d.(d_ch).ModuleType1.R0_matCell;
        % end

    else
        d.(d_ch).ModuleType1.R0_matCell = x_RC.("R0_mat_"+d_ch);
        d.(d_ch).R_loss = R_loss;
    end

    %% ModuleType1_original
    
    %Stock Simulink Values:
    %ModuleType1.SOC_vec = [0, .1, .25, .5, .75, .9, 1]; % Vector of state-of-charge values, SOC
    %ModuleType1.T_vec = [278, 293, 313]; % Vector of temperatures, T, K
    %ModuleType1.V0_mat = [3.49, 3.5, 3.51; 3.55, 3.57, 3.56; 3.62, 3.63, 3.64; 3.71, 3.71, 3.72; 3.91, 3.93, 3.94; 4.07, 4.08, 4.08; 4.19, 4.19, 4.19]; % Open-circuit voltage, V0(SOC,T), V
    
    %ModuleType1.R0_mat = [.0117, .0085, .009; .011, .0085, .009; .0114, .0087, .0092; .0107, .0082, .0088; .0107, .0083, .0091; .0113, .0085, .0089; .0116, .0085, .0089]; % Terminal resistance, R0(SOC,T), Ohm
    
    d.(d_ch).ModuleType1.AHCell = battery_GD.SysInfo.C_sys; % Cell capacity, AH, A*hr
    
    %ModuleType1.R1_mat = [.0109, .0029, .0013; .0069, .0024, .0012; .0047, .0026, .0013; .0034, .0016, .001; .0033, .0023, .0014; .0033, .0018, .0011; .0028, .0017, .0011]; % First polarization resistance, R1(SOC,T), Ohm
    %ModuleType1.tau1_mat = [20, 36, 39; 31, 45, 39; 109, 105, 61; 36, 29, 26; 59, 77, 67; 40, 33, 29; 25, 39, 33]; % First time constant, tau1(SOC,T), s
    %ModuleType1.R2_mat = [.0109, .0029, .0013; .0069, .0024, .0012; .0047, .0026, .0013; .0034, .0016, .001; .0033, .0023, .0014; .0033, .0018, .0011; .0028, .0017, .0011]; % Second polarization resistance, R2(SOC,T), Ohm
    %ModuleType1.tau2_mat = [20, 36, 39; 31, 45, 39; 109, 105, 61; 36, 29, 26; 59, 77, 67; 40, 33, 29; 25, 39, 33]; % Second time constant, tau2(SOC,T), s

    % Stock Simscape values to keep
    d.(d_ch).ModuleType1.thermal_massCell = 100; % Thermal mass, J/K
    d.(d_ch).ModuleType1.CoolantResistance = 1.2; % Cell level coolant thermal path resistance, K/W
    d.(d_ch).ModuleType1.AmbientResistance = 25; % Cell level ambient thermal path resistance, K/W

    d.(d_ch).ModuleType1.V_rangeCell = [0, inf]; % Terminal voltage operating range [Min Max], V
    
    d.(d_ch).ModuleType1.CellBalancingClosedResistance = 0.01; % Cell balancing switch closed resistance, Ohm
    d.(d_ch).ModuleType1.CellBalancingOpenConductance = 1e-8; % Cell balancing switch open conductance, 1/Ohm
    d.(d_ch).ModuleType1.CellBalancingThreshold = 0.5; % Cell balancing switch operation threshold
    d.(d_ch).ModuleType1.CellBalancingResistance = 50; % Cell balancing shunt resistance, Ohm
    
    disp("battery_initializer: Successfully initialized ModuleType1")
    
    %% ParallelAssemblyType1 ------ Test how Parallel Assemblence influences these Values
    %{
    ParallelAssemblyType1.SOC_vec = [0, .1, .25, .5, .75, .9, 1]; % Vector of state-of-charge values, SOC
    ParallelAssemblyType1.T_vec = [278, 293, 313]; % Vector of temperatures, T, K
    ParallelAssemblyType1.V0_mat = [3.49, 3.5, 3.51; 3.55, 3.57, 3.56; 3.62, 3.63, 3.64; 3.71, 3.71, 3.72; 3.91, 3.93, 3.94; 4.07, 4.08, 4.08; 4.19, 4.19, 4.19]; % Open-circuit voltage, V0(SOC,T), V
    ParallelAssemblyType1.V_range = [0, inf]; % Terminal voltage operating range [Min Max], V
    ParallelAssemblyType1.R0_mat = [.0117, .0085, .009; .011, .0085, .009; .0114, .0087, .0092; .0107, .0082, .0088; .0107, .0083, .0091; .0113, .0085, .0089; .0116, .0085, .0089]; % Terminal resistance, R0(SOC,T), Ohm
    ParallelAssemblyType1.AH = 27; % Cell capacity, AH, A*hr
    ParallelAssemblyType1.R1_mat = [.0109, .0029, .0013; .0069, .0024, .0012; .0047, .0026, .0013; .0034, .0016, .001; .0033, .0023, .0014; .0033, .0018, .0011; .0028, .0017, .0011]; % First polarization resistance, R1(SOC,T), Ohm
    ParallelAssemblyType1.tau1_mat = [20, 36, 39; 31, 45, 39; 109, 105, 61; 36, 29, 26; 59, 77, 67; 40, 33, 29; 25, 39, 33]; % First time constant, tau1(SOC,T), s
    ParallelAssemblyType1.R2_mat = [.0109, .0029, .0013; .0069, .0024, .0012; .0047, .0026, .0013; .0034, .0016, .001; .0033, .0023, .0014; .0033, .0018, .0011; .0028, .0017, .0011]; % Second polarization resistance, R2(SOC,T), Ohm
    ParallelAssemblyType1.tau2_mat = [20, 36, 39; 31, 45, 39; 109, 105, 61; 36, 29, 26; 59, 77, 67; 40, 33, 29; 25, 39, 33]; % Second time constant, tau2(SOC,T), s
    ParallelAssemblyType1.thermal_mass = 100; % Thermal mass, J/K
    ParallelAssemblyType1.CoolantResistance = 1.2; % Cell level coolant thermal path resistance, K/W
    ParallelAssemblyType1.AmbientResistance = 25; % Cell level ambient thermal path resistance, K/W
    ParallelAssemblyType1.CellBalancingClosedResistance = 0.01; % Cell balancing switch closed resistance, Ohm
    ParallelAssemblyType1.CellBalancingOpenConductance = 1e-8; % Cell balancing switch open conductance, 1/Ohm
    ParallelAssemblyType1.CellBalancingThreshold = 0.5; % Cell balancing switch operation threshold
    ParallelAssemblyType1.CellBalancingResistance = 50; % Cell balancing shunt resistance, Ohm
    %}
    
    %Preliminary Soultion:
    for i = h.ParallelAssemblyType1_fieldnames
        d.(d_ch).ParallelAssemblyType1.(i) = d.(d_ch).ModuleType1.(i);
    end
    
    disp("battery_initializer: Successfully initialized ParallelAssemblyType1")
    
    %% ModuleAssembly1.Module1
    h.socCellModule = repmat(config.simulation_SOC_start,battery_GD.SysInfo.num_serial_mods_sys,1);
    h.socParallelAssembly = h.socCellModule;
    if config.simulation_SOC_deviation > 0
        if config.simulation_SOC_start >= 1 || config.simulation_SOC_start <= 0
            error("battery_initializer: Can't deviate with the chosen SOC")
        else
            h.socCellModule = h.socCellModule * (1 + randn(battery_GD.SysInfo.num_serial_mods_sys,1) * config.simulation_SOC_deviation);
            h.socParallelAssembly = mean(h.socCellModule);
        end
    end
    
    if config.simulation_cell_grouping == "Lumped"
        h.n_c_m = 1;
    else
        h.n_c_m = battery_GD.ModInfo.num_cells_mod;
    end
    
    h.n_s_c_m = battery_GD.ModInfo.num_serial_cells_mod;
    h.n_s_m_s = battery_GD.SysInfo.num_serial_mods_sys;
    
    for j = 1:h.n_s_m_s
        if h.n_s_m_s > 10 && j < 10
            h.modules = sprintf("Module0%d", j);
        else
            h.modules = sprintf("Module%d", j);
        end

        for i = 1:battery_GD.SysInfo.num_serial_mods_sys
            d.(d_ch).("ModuleAssembly"+string(i)).(h.modules).iCell = repmat(0, h.n_c_m, 1); % Cell model current (positive in), A        %standard value
            d.(d_ch).("ModuleAssembly"+string(i)).(h.modules).vCell = repmat(0, h.n_c_m, 1); % Cell model terminal voltage, V             %standard value
            d.(d_ch).("ModuleAssembly"+string(i)).(h.modules).socCell = repmat(h.socCellModule(i), h.n_c_m, 1); % Cell model state of charge    
            d.(d_ch).("ModuleAssembly"+string(i)).(h.modules).numCyclesCell = repmat(0, h.n_c_m, 1); % Cell model discharge cycles        %standard value
            d.(d_ch).("ModuleAssembly"+string(i)).(h.modules).temperatureCell = repmat(config.vehicle.ambient_temp, h.n_c_m, 1); % Cell model temperature, K   %standard value
            d.(d_ch).("ModuleAssembly"+string(i)).(h.modules).vParallelAssembly = repmat(0, h.n_s_c_m, 1); % Parallel Assembly Voltage, V        %standard value
            d.(d_ch).("ModuleAssembly"+string(i)).(h.modules).socParallelAssembly = repmat(h.socCellModule(i), h.n_s_c_m, 1); % Parallel Assembly state of charge
        end
    end
    disp("battery_initializer: Successfully initialized all ModuleAssemblys")
end

clear d_ch
%% save Initialization
ch = d.ch;
save(pathes.bat.Batterymodel+"/+"+config.libraryname+"/"+config.libraryname+"_param_ch.mat","-struct","ch") % save Charge Initilaization

dch = d.dch;
save(pathes.bat.Batterymodel+"/+"+config.libraryname+"/"+config.libraryname+"_param_dch.mat","-struct","dch") % save Discharge Initilaization

disp("battery_initializer: Successfully saved Simscape Library initialization")
clear ch dch

%% Other Preperations for the simualtion
% Surface for natural heat convection (worst case apporach excluding "_BTMS")
battery_GD.SysInfo.dim_surface = ...
        2* (battery_GD.SysInfo.dim_x_sys * battery_GD.SysInfo.dim_y_sys +...
        battery_GD.SysInfo.dim_x_sys * battery_GD.SysInfo.dim_z_sys +...
        battery_GD.SysInfo.dim_y_sys * battery_GD.SysInfo.dim_z_sys);

% Voltage at simulation start
h.x = 1:size(d.dch.ModuleType1.V0_matCell,1);

h.y = d.dch.ModuleType1.V0_matCell(:,ceil(size(x_RC.T_vec,2)));
h.xq = 1 + (size(d.dch.ModuleType1.V0_matCell,1)-1) * config.simulation_SOC_start;
battery_GD.SysInfo.U0_start = battery_GD.ModInfo.num_serial_cells_mod * ...
        battery_GD.SysInfo.num_serial_mods_sys * ...
        interp1(h.x, h.y , h.xq , 'linear');


%% assign values into workspace

%{
notassign2workspace = ["battery_GD", "config", "pathes", "notassign2workspace"]; %remove help variables
assign2workspace = setdiff(who,notassign2workspace);

for i=1:size(assign2workspace,1)
   %assignin("caller",assign2workspace(i),eval(assign2workspace(i)));
   batteryInitialization.(assign2workspace(i)) = eval(assign2workspace(i));
end
%}

% Suppress MATLAB editor message regarding readability of repmat
%#ok<*REPMAT>

end
