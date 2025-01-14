function [battery_LossProbe] = battery_LossProbeAdder(battery_GD)
%% battery_LossProbeAdder v1.3
% Includes:
% - All loss components (Cell, Balancing Resistor, Balancing Switch, Coolant)
% - Scaling by system and module parameters
% - Two output ports positioned below sink_bat_loss_out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Initialization
battery_LossProbe = 0; % Flag of function state
model = evalin('base', 'model');
subsys = [model '/Variant Subsystem Battery/Battery Global Drive/Battery'];

% Block names
ProbeBlockName = [subsys '/Pack1/ModuleAssembly1/LossProbe'];
BoundBlockName = [subsys '/Pack1/ModuleAssembly1/Module01'];
SinkBlockName = [subsys '/Pack1/ModuleAssembly1/sink_bat_loss_out'];
SumBlockName = [subsys '/Pack1/ModuleAssembly1/Sum_Losses'];
GainBlockName = [subsys '/Pack1/ModuleAssembly1/Scale_TotalLoss'];
OutportBlockName = [subsys '/Pack1/ModuleAssembly1/Out_bat_loss'];
OutportBlockName2 = [subsys '/Pack1/Out_bat_loss'];
SysDepthName = [subsys '/Pack1/ModuleAssembly1'];

% Block positions
ProbeBlockPosition = [-200, 100, -50, 250];      % LossProbe enlarged and moved left
SumBlockPosition = [-40, 140, 0, 170];          % Sum_Losses moved left
GainBlockPosition = [40, 140, 80, 170];         % Scale_TotalLoss moved left
SinkBlockPosition = [120, 140, 200, 170];       % sink_bat_loss_out
OutportBlockPosition = [140, 200, 180, 220];    % Out_bat_loss positioned below sink_bat_loss_out
NewOutportBlockPosition = [140, 300, 180, 320]; % New_Out_bat_loss positioned below Out_bat_loss

load_system(model)

%% Clear Old Blocks
try delete_line(SysDepthName, 'LossProbe/1', 'Sum_Losses/1'); catch, end
try delete_line(SysDepthName, 'Sum_Losses/1', 'Scale_TotalLoss/1'); catch, end
try delete_line(SysDepthName, 'Scale_TotalLoss/1', 'sink_bat_loss_out/1'); catch, end
try delete_line(SysDepthName, 'sink_bat_loss_out/1', 'Out_bat_loss/1'); catch, end
try delete_line(SysDepthName, 'Sum_Losses/1', 'New_Out_bat_loss/1'); catch, end
try delete_block(ProbeBlockName); catch, end
try delete_block(SumBlockName); catch, end
try delete_block(GainBlockName); catch, end
try delete_block(SinkBlockName); catch, end
try delete_block(OutportBlockName); catch, end
try delete_block(OutportBlockName2); catch, end

%% Set up new Probing
% Add Probe block
probeBlock = add_block('nesl_utility/Probe', ProbeBlockName);
simscape.probe.setBoundBlock(ProbeBlockName, BoundBlockName);
simscape.probe.setVariables(ProbeBlockName, ...
    ["Cell1.power_dissipated", ...
     "CoolantResistor.Q", ...
     "balancingResistor.power_dissipated", ...
     "balancingSwitch.power_dissipated"]);
set_param(ProbeBlockName, 'Position', ProbeBlockPosition);

%% Add Sum Block
sumBlock = add_block('simulink/Math Operations/Sum', SumBlockName, ...
                     'Inputs', '++++', 'Position', SumBlockPosition);

%% Add Gain Block
% Define scaling factors
num_parallel_mods_sys = battery_GD.SysInfo.num_parallel_mods_sys;
num_serial_mods_sys = battery_GD.SysInfo.num_serial_mods_sys;
num_parallel_cells_mod = battery_GD.ModInfo.num_parallel_cells_mod;
num_serial_cells_mod = battery_GD.ModInfo.num_serial_cells_mod;

total_gain = num_parallel_mods_sys * num_serial_mods_sys * ...
             num_parallel_cells_mod * num_serial_cells_mod;

gainBlock = add_block('simulink/Math Operations/Gain', GainBlockName, ...
                      'Gain', num2str(total_gain), ...
                      'Position', GainBlockPosition);

%% Add To Workspace Block
toworkspaceBlock = add_block('simulink/Sinks/To Workspace', SinkBlockName, ...
                             'VariableName', 'bat_loss_out', ...
                             'Position', SinkBlockPosition);

%% Add Outport Blocks
outportBlock = add_block('simulink/Ports & Subsystems/Out1', OutportBlockName, ...
                         'Position', OutportBlockPosition);

newOutportBlock = add_block('simulink/Ports & Subsystems/Out1', OutportBlockName2, ...
                            'Position', NewOutportBlockPosition);

%% Connect Blocks
% Connect Probe to Sum block
add_line(SysDepthName, 'LossProbe/1', 'Sum_Losses/1'); % Cell1.power_dissipated
add_line(SysDepthName, 'LossProbe/2', 'Sum_Losses/2'); % CoolantResistor.Q
add_line(SysDepthName, 'LossProbe/3', 'Sum_Losses/3'); % BalancingResistor
add_line(SysDepthName, 'LossProbe/4', 'Sum_Losses/4'); % BalancingSwitch

% Connect Sum block to Gain block
add_line(SysDepthName, 'Sum_Losses/1', 'Scale_TotalLoss/1');

% Connect Gain block to To Workspace block
add_line(SysDepthName, 'Scale_TotalLoss/1', 'sink_bat_loss_out/1');

% Connect Gain block to Outport block
add_line(SysDepthName, 'Scale_TotalLoss/1', 'Out_bat_loss/1');

% Connect Sum block to New Outport block
add_line([subsys '/Pack1'], 'ModuleAssembly1/8', 'Out_bat_loss/1');

%Connect Pack Output to Outport block
add_line(subsys, 'Pack1/8', 'battery_loss/1');

%% Save and Close System
save_system(model)
close_system(model)
battery_LossProbe = 1; % Flag of function state   
fprintf("Successfully placed Probe for Battery Loss with scaling and added two Outport blocks\n")

end
