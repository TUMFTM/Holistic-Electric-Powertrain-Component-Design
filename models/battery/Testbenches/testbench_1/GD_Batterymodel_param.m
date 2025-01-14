%% Battery parameters

%% ModuleType1
ModuleType1.SOC_vec = [0, .1, .25, .5, .75, .9, 1]; % Vector of state-of-charge values, SOC
ModuleType1.V0_vec = [3.5057, 3.566, 3.6337, 3.7127, 3.9259, 4.0777, 4.1928]; % Open-circuit voltage, V0(SOC), V
ModuleType1.V_range = [0, inf]; % Terminal voltage operating range [Min Max], V
ModuleType1.R0_vec = [.0085, .0085, .0087, .0082, .0083, .0085, .0085]; % Terminal resistance, R0(SOC), Ohm
ModuleType1.AH = 27; % Cell capacity, AH, A*hr

%% ParallelAssemblyType1
ParallelAssemblyType1.SOC_vec = [0, .1, .25, .5, .75, .9, 1]; % Vector of state-of-charge values, SOC
ParallelAssemblyType1.V0_vec = [3.5057, 3.566, 3.6337, 3.7127, 3.9259, 4.0777, 4.1928]; % Open-circuit voltage, V0(SOC), V
ParallelAssemblyType1.V_range = [0, inf]; % Terminal voltage operating range [Min Max], V
ParallelAssemblyType1.R0_vec = [.0085, .0085, .0087, .0082, .0083, .0085, .0085]; % Terminal resistance, R0(SOC), Ohm
ParallelAssemblyType1.AH = 27; % Cell capacity, AH, A*hr

%% Battery initial targets

%% ModuleAssembly1.Module1
ModuleAssembly1.Module1.iCellModel = 0; % Cell model current (positive in), A
ModuleAssembly1.Module1.vCellModel = 0; % Cell model terminal voltage, V
ModuleAssembly1.Module1.socCellModel = 1; % Cell model state of charge
ModuleAssembly1.Module1.numCyclesCellModel = 0; % Cell model discharge cycles
ModuleAssembly1.Module1.vParallelAssembly = repmat(0, 14, 1); % Parallel Assembly Voltage, V
ModuleAssembly1.Module1.socParallelAssembly = repmat(1, 14, 1); % Parallel Assembly state of charge

%% ModuleAssembly2.Module1
ModuleAssembly2.Module1.iCellModel = 0; % Cell model current (positive in), A
ModuleAssembly2.Module1.vCellModel = 0; % Cell model terminal voltage, V
ModuleAssembly2.Module1.socCellModel = 0.8; % Cell model state of charge
ModuleAssembly2.Module1.numCyclesCellModel = 0; % Cell model discharge cycles
ModuleAssembly2.Module1.vParallelAssembly = repmat(0, 14, 1); % Parallel Assembly Voltage, V
ModuleAssembly2.Module1.socParallelAssembly = repmat(1, 14, 1); % Parallel Assembly state of charge

% Suppress MATLAB editor message regarding readability of repmat
%#ok<*REPMAT>
