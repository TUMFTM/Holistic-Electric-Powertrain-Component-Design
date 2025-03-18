%% Battery parameters

%% ModuleType1
ModuleType1.SOC_vecCell = [0, .1, .25, .5, .75, .9, 1]; % Vector of state-of-charge values, SOC
ModuleType1.T_vecCell = [278, 293, 313]; % Vector of temperatures, T, K
ModuleType1.V0_matCell = [3.49, 3.5, 3.51; 3.55, 3.57, 3.56; 3.62, 3.63, 3.64; 3.71, 3.71, 3.72; 3.91, 3.93, 3.94; 4.07, 4.08, 4.08; 4.19, 4.19, 4.19]; % Open-circuit voltage, V0(SOC,T), V
ModuleType1.V_rangeCell = [0, inf]; % Terminal voltage operating range [Min Max], V
ModuleType1.R0_matCell = [.0117, .0085, .009; .011, .0085, .009; .0114, .0087, .0092; .0107, .0082, .0088; .0107, .0083, .0091; .0113, .0085, .0089; .0116, .0085, .0089]; % Terminal resistance, R0(SOC,T), Ohm
ModuleType1.AHCell = 27; % Cell capacity, AH, A*hr
ModuleType1.R1_matCell = [.0109, .0029, .0013; .0069, .0024, .0012; .0047, .0026, .0013; .0034, .0016, .001; .0033, .0023, .0014; .0033, .0018, .0011; .0028, .0017, .0011]; % First polarization resistance, R1(SOC,T), Ohm
ModuleType1.tau1_matCell = [20, 36, 39; 31, 45, 39; 109, 105, 61; 36, 29, 26; 59, 77, 67; 40, 33, 29; 25, 39, 33]; % First time constant, tau1(SOC,T), s
ModuleType1.R2_matCell = [.0109, .0029, .0013; .0069, .0024, .0012; .0047, .0026, .0013; .0034, .0016, .001; .0033, .0023, .0014; .0033, .0018, .0011; .0028, .0017, .0011]; % Second polarization resistance, R2(SOC,T), Ohm
ModuleType1.tau2_matCell = [20, 36, 39; 31, 45, 39; 109, 105, 61; 36, 29, 26; 59, 77, 67; 40, 33, 29; 25, 39, 33]; % Second time constant, tau2(SOC,T), s
ModuleType1.thermal_massCell = 100; % Thermal mass, J/K
ModuleType1.CoolantResistance = 1.2; % Cell level coolant thermal path resistance, K/W
ModuleType1.AmbientResistance = 25; % Cell level ambient thermal path resistance, K/W
ModuleType1.CellBalancingClosedResistance = 0.01; % Cell balancing switch closed resistance, Ohm
ModuleType1.CellBalancingOpenConductance = 1e-8; % Cell balancing switch open conductance, 1/Ohm
ModuleType1.CellBalancingThreshold = 0.5; % Cell balancing switch operation threshold
ModuleType1.CellBalancingResistance = 50; % Cell balancing shunt resistance, Ohm

%% ParallelAssemblyType1
ParallelAssemblyType1.SOC_vecCell = [0, .1, .25, .5, .75, .9, 1]; % Vector of state-of-charge values, SOC
ParallelAssemblyType1.T_vecCell = [278, 293, 313]; % Vector of temperatures, T, K
ParallelAssemblyType1.V0_matCell = [3.49, 3.5, 3.51; 3.55, 3.57, 3.56; 3.62, 3.63, 3.64; 3.71, 3.71, 3.72; 3.91, 3.93, 3.94; 4.07, 4.08, 4.08; 4.19, 4.19, 4.19]; % Open-circuit voltage, V0(SOC,T), V
ParallelAssemblyType1.V_rangeCell = [0, inf]; % Terminal voltage operating range [Min Max], V
ParallelAssemblyType1.R0_matCell = [.0117, .0085, .009; .011, .0085, .009; .0114, .0087, .0092; .0107, .0082, .0088; .0107, .0083, .0091; .0113, .0085, .0089; .0116, .0085, .0089]; % Terminal resistance, R0(SOC,T), Ohm
ParallelAssemblyType1.AHCell = 27; % Cell capacity, AH, A*hr
ParallelAssemblyType1.R1_matCell = [.0109, .0029, .0013; .0069, .0024, .0012; .0047, .0026, .0013; .0034, .0016, .001; .0033, .0023, .0014; .0033, .0018, .0011; .0028, .0017, .0011]; % First polarization resistance, R1(SOC,T), Ohm
ParallelAssemblyType1.tau1_matCell = [20, 36, 39; 31, 45, 39; 109, 105, 61; 36, 29, 26; 59, 77, 67; 40, 33, 29; 25, 39, 33]; % First time constant, tau1(SOC,T), s
ParallelAssemblyType1.R2_matCell = [.0109, .0029, .0013; .0069, .0024, .0012; .0047, .0026, .0013; .0034, .0016, .001; .0033, .0023, .0014; .0033, .0018, .0011; .0028, .0017, .0011]; % Second polarization resistance, R2(SOC,T), Ohm
ParallelAssemblyType1.tau2_matCell = [20, 36, 39; 31, 45, 39; 109, 105, 61; 36, 29, 26; 59, 77, 67; 40, 33, 29; 25, 39, 33]; % Second time constant, tau2(SOC,T), s
ParallelAssemblyType1.thermal_massCell = 100; % Thermal mass, J/K
ParallelAssemblyType1.CoolantResistance = 1.2; % Cell level coolant thermal path resistance, K/W
ParallelAssemblyType1.AmbientResistance = 25; % Cell level ambient thermal path resistance, K/W
ParallelAssemblyType1.CellBalancingClosedResistance = 0.01; % Cell balancing switch closed resistance, Ohm
ParallelAssemblyType1.CellBalancingOpenConductance = 1e-8; % Cell balancing switch open conductance, 1/Ohm
ParallelAssemblyType1.CellBalancingThreshold = 0.5; % Cell balancing switch operation threshold
ParallelAssemblyType1.CellBalancingResistance = 50; % Cell balancing shunt resistance, Ohm

%% Battery initial targets

%% ModuleAssembly1.Module1
ModuleAssembly1.Module1.iCell = 0; % Cell current (positive in), A
ModuleAssembly1.Module1.vCell = 0; % Cell terminal voltage, V
ModuleAssembly1.Module1.socCell = 1; % Cell state of charge
ModuleAssembly1.Module1.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly1.Module1.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly1.Module1.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly1.Module1.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly1.Module2
ModuleAssembly1.Module2.iCell = 0; % Cell current (positive in), A
ModuleAssembly1.Module2.vCell = 0; % Cell terminal voltage, V
ModuleAssembly1.Module2.socCell = 1; % Cell state of charge
ModuleAssembly1.Module2.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly1.Module2.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly1.Module2.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly1.Module2.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly1.Module3
ModuleAssembly1.Module3.iCell = 0; % Cell current (positive in), A
ModuleAssembly1.Module3.vCell = 0; % Cell terminal voltage, V
ModuleAssembly1.Module3.socCell = 1; % Cell state of charge
ModuleAssembly1.Module3.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly1.Module3.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly1.Module3.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly1.Module3.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly1.Module4
ModuleAssembly1.Module4.iCell = 0; % Cell current (positive in), A
ModuleAssembly1.Module4.vCell = 0; % Cell terminal voltage, V
ModuleAssembly1.Module4.socCell = 1; % Cell state of charge
ModuleAssembly1.Module4.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly1.Module4.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly1.Module4.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly1.Module4.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly1.Module5
ModuleAssembly1.Module5.iCell = 0; % Cell current (positive in), A
ModuleAssembly1.Module5.vCell = 0; % Cell terminal voltage, V
ModuleAssembly1.Module5.socCell = 1; % Cell state of charge
ModuleAssembly1.Module5.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly1.Module5.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly1.Module5.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly1.Module5.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly1.Module6
ModuleAssembly1.Module6.iCell = 0; % Cell current (positive in), A
ModuleAssembly1.Module6.vCell = 0; % Cell terminal voltage, V
ModuleAssembly1.Module6.socCell = 1; % Cell state of charge
ModuleAssembly1.Module6.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly1.Module6.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly1.Module6.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly1.Module6.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly2.Module1
ModuleAssembly2.Module1.iCell = 0; % Cell current (positive in), A
ModuleAssembly2.Module1.vCell = 0; % Cell terminal voltage, V
ModuleAssembly2.Module1.socCell = 1; % Cell state of charge
ModuleAssembly2.Module1.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly2.Module1.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly2.Module1.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly2.Module1.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly2.Module2
ModuleAssembly2.Module2.iCell = 0; % Cell current (positive in), A
ModuleAssembly2.Module2.vCell = 0; % Cell terminal voltage, V
ModuleAssembly2.Module2.socCell = 1; % Cell state of charge
ModuleAssembly2.Module2.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly2.Module2.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly2.Module2.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly2.Module2.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly2.Module3
ModuleAssembly2.Module3.iCell = 0; % Cell current (positive in), A
ModuleAssembly2.Module3.vCell = 0; % Cell terminal voltage, V
ModuleAssembly2.Module3.socCell = 1; % Cell state of charge
ModuleAssembly2.Module3.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly2.Module3.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly2.Module3.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly2.Module3.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly2.Module4
ModuleAssembly2.Module4.iCell = 0; % Cell current (positive in), A
ModuleAssembly2.Module4.vCell = 0; % Cell terminal voltage, V
ModuleAssembly2.Module4.socCell = 1; % Cell state of charge
ModuleAssembly2.Module4.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly2.Module4.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly2.Module4.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly2.Module4.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly2.Module5
ModuleAssembly2.Module5.iCell = 0; % Cell current (positive in), A
ModuleAssembly2.Module5.vCell = 0; % Cell terminal voltage, V
ModuleAssembly2.Module5.socCell = 1; % Cell state of charge
ModuleAssembly2.Module5.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly2.Module5.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly2.Module5.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly2.Module5.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

%% ModuleAssembly2.Module6
ModuleAssembly2.Module6.iCell = 0; % Cell current (positive in), A
ModuleAssembly2.Module6.vCell = 0; % Cell terminal voltage, V
ModuleAssembly2.Module6.socCell = 1; % Cell state of charge
ModuleAssembly2.Module6.numCyclesCell = 0; % Cell discharge cycles
ModuleAssembly2.Module6.temperatureCell = 298.15; % Cell temperature, K
ModuleAssembly2.Module6.vParallelAssembly = repmat(0, 18, 1); % Parallel Assembly Voltage, V
ModuleAssembly2.Module6.socParallelAssembly = repmat(1, 18, 1); % Parallel Assembly state of charge

% Suppress MATLAB editor message regarding readability of repmat
%#ok<*REPMAT>
