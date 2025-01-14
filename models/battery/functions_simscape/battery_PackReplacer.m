function [battery_PackReplacement] = battery_PackReplacer()
%% battery_PackReplacer v0.1
% This function was created as part of Global Drive 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                           Versions                                    %
% v0.1 - initial setup                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        known Problems                                 %
% - none                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation
battery_PackReplacement = 0; %Flag of function state
model = evalin('base', 'model');
% model = 'BEV_Laengsdynamikmodell_v2';
subsys = [model '/Variant Subsystem Battery/Battery Global Drive/Battery'];
% model = '03_Messdaten/WLTP/Nachsimulation/BEV_Laengsdynamikmodell_ID3_Test';
% model = 'Efficiency_Testbench'; % Test/Debug model

lib.Path = 'models/battery/Batterymodel/GD_Batterymodel_NSGAII_1';
lib.Name = 'GD_Batterymodel_NSGAII_1';

load_system(model)
load_system(lib.Path)

% replace_block(model, Pack.Name, [lib.Name '/Pack1']); % Does not work - invalid block type

blockList = find_system([model '/Variant Subsystem Battery/Battery Global Drive/Battery'], 'SearchDepth', 1, 'Type', 'Block');
Pack.Name = [model '/Variant Subsystem Battery/Battery Global Drive/Battery/Pack1'];
Pack.LineHandle = get_param(Pack.Name,"LineHandles");
Pack.PortHandle = get_param(Pack.Name,"PortHandles");
Pack.Position = get_param(Pack.Name,"Position");
Pack.OutportName = get_param(Pack.PortHandle.Outport,'Name');

% Get remaining block Names
for i = 1:numel(blockList)
    % oben
    if ~isempty(strfind(blockList{i}, "Controlled Current"))
        CCS.Name = blockList{i};
    % rechts
    elseif ~isempty(strfind(blockList{i}, "Mean_i"))
        Mean_i.Name = blockList{i};
    elseif ~isempty(strfind(blockList{i}, 'Terminator_numcyc'))
        Term_numcyc.Name = blockList{i};
    elseif ~isempty(strfind(blockList{i}, "Terminator_socC"))
        Term_socC.Name = blockList{i};
    elseif ~isempty(strfind(blockList{i}, "Mean_SOC"))
        Mean_SOC.Name = blockList{i};
    elseif ~isempty(strfind(blockList{i}, "Mean_temp"))
        Mean_temp.Name = blockList{i};
    elseif ~isempty(strfind(blockList{i}, "Mean_v"))
        Mean_v.Name = blockList{i};
    elseif ~isempty(strfind(blockList{i}, "Passive"))
        Passive.Name = blockList{i};
    % unten
    elseif ~isempty(strfind(blockList{i}, "Solver"))
        Solver.Name = blockList{i};
    elseif ~isempty(strfind(blockList{i}, "Convective"))
        Convective.Name = blockList{i};
    elseif ~isempty(strfind(blockList{i}, "Parallel"))
        ParallelChannels.Name = blockList{i};
    % links
    elseif ~isempty(strfind(blockList{i}, "Delay_PB"))
        UnitDelay.Name = blockList{i};
    end
end

% oben
% Controlled Current Source
CCS.LineHandle = get_param(CCS.Name,"LineHandles");
CCS.PortHandle = get_param(CCS.Name,"PortHandles");

% rechts
% Mean_i
Mean_i.PortHandle = get_param(Mean_i.Name,"PortHandles");
Mean_i.InportName = get_param(Mean_i.PortHandle.Inport,'Name');
% Term_numcyc
Term_numcyc.PortHandle = get_param(Term_numcyc.Name,"PortHandles");
Term_numcyc.InportName = get_param(Term_numcyc.PortHandle.Inport,'Name');
% Term_socC
Term_socC.PortHandle = get_param(Term_socC.Name,"PortHandles");
Term_socC.InportName = get_param(Term_socC.PortHandle.Inport,'Name');
% Mean_SOC
Mean_SOC.PortHandle = get_param(Mean_SOC.Name,"PortHandles");
Mean_SOC.InportName = get_param(Mean_SOC.PortHandle.Inport,'Name');
% Mean_temp
Mean_temp.PortHandle = get_param(Mean_temp.Name,"PortHandles");
Mean_temp.InportName = get_param(Mean_temp.PortHandle.Inport,'Name');
% Mean_v
Mean_v.PortHandle = get_param(Mean_v.Name,"PortHandles");
Mean_v.InportName = get_param(Mean_v.PortHandle.Inport,'Name');
%Passive
Passive.PortHandle = get_param(Passive.Name,"PortHandles");
Passive.InportName = get_param(Passive.PortHandle.Inport,'Name');

% unten
% SolverConfiguration
Solver.LineHandle = get_param(Solver.Name,"LineHandles");
Solver.PortHandle = get_param(Solver.Name,"PortHandles");

% Convective
Convective.LineHandle = get_param(Convective.Name,"LineHandles");
Convective.PortHandle = get_param(Convective.Name,"PortHandles");

% ParallelChannels
ParallelChannels.LineHandle = get_param(ParallelChannels.Name,"LineHandles");
ParallelChannels.PortHandle = get_param(ParallelChannels.Name,"PortHandles");
%

% links
% Unity Delay_PB
UnitDelay.LineHandle = get_param(UnitDelay.Name,"LineHandles");
UnitDelay.PortHandle = get_param(UnitDelay.Name,"PortHandles");
%

% if battery_GD.SysInfo.num_serial_mods_sys < 11
%     BoundBlockName = [model '/Pack1/ModuleAssembly1/Module1'];
% else
%     BoundBlockName = [model '/Pack1/ModuleAssembly1/Module01'];
% end
% SinkBlockName = [model '/Pack1/ModuleAssembly1/sink_bat_loss_out'];
% SysDepthName = [model '/Pack1'];
% SinkBlockPostion = [165   140   255   170]; 
% % % outport_MA1Name = 'Efficiency_Testbench/Pack1/ModuleAssembly1/LossOutport_MA';
% % % outport_PackName = 'Efficiency_Testbench/Pack1/LossOutport_Pack';

load_system(model)
load_system(lib.Path)

%% Clear Old junk, to prevent crashes

for i = 1:13
    try
        switch i
        %oben
            case 1    
                delete_line(CCS.LineHandle.LConn); % Works
        % rechts    
            case 2
                delete_line(Pack.LineHandle.Outport(1)) % Works
            case 3
                delete_line(Pack.LineHandle.Outport(2)) % Works
            case 4
                delete_line(Pack.LineHandle.Outport(3)) % Works
            case 5
                delete_line(Pack.LineHandle.Outport(4)) % Works
            case 6
                delete_line(Pack.LineHandle.Outport(5)) % Works
            case 7
                delete_line(Pack.LineHandle.Outport(6)) % Works
            case 8
                delete_line(Pack.LineHandle.Outport(7)) % Works
            case 9
                delete_line(Pack.LineHandle.Outport(8)) % Works
        %unten
            case 10
                delete_line(Pack.LineHandle.RConn(1)); % Works
            case 11
                delete_line(Pack.LineHandle.RConn(2)); % Works
            case 12
                delete_line(Pack.LineHandle.RConn(3)); % Works
        %links
            case 13
                delete_line(Pack.LineHandle.Inport) % Works
        end
    catch
    end
end

try
    delete_block(Pack.Name); % Works
catch
end

%% Set up new Pack

add_block([lib.Name '/Pack1'], Pack.Name); % Works
close_system(lib.Name)
set_param(Pack.Name,'Position', Pack.Position); % Works
% Pack.LineHandle = get_param(Pack.Name,"LineHandles");
Pack.PortHandle = get_param(Pack.Name,"PortHandles");

% Add Lines
% oben
add_line(subsys, CCS.PortHandle.LConn,Pack.PortHandle.LConn, 'autorouting', 'on') % Works

% rechts
add_line(subsys, Pack.PortHandle.Outport(1),Mean_i.PortHandle.Inport, 'autorouting', 'on') % Works

add_line(subsys, Pack.PortHandle.Outport(2),Term_numcyc.PortHandle.Inport, 'autorouting', 'on') % Works
add_line(subsys, Pack.PortHandle.Outport(3),Term_socC.PortHandle.Inport, 'autorouting', 'on') % Works
add_line(subsys, Pack.PortHandle.Outport(4),Mean_SOC.PortHandle.Inport, 'autorouting', 'on') % Works
add_line(subsys, Pack.PortHandle.Outport(5),Mean_temp.PortHandle.Inport, 'autorouting', 'on') % Works
add_line(subsys, Pack.PortHandle.Outport(6),Mean_v.PortHandle.Inport, 'autorouting', 'on') % Works
add_line(subsys, Pack.PortHandle.Outport(7),Passive.PortHandle.Inport(1), 'autorouting', 'on') % Works

% unten
add_line(subsys, Pack.PortHandle.RConn(1),Solver.PortHandle.RConn, 'autorouting', 'on') % Works
add_line(subsys, Pack.PortHandle.RConn(2),Convective.PortHandle.LConn, 'autorouting', 'on') % Works
add_line(subsys, Pack.PortHandle.RConn(3),ParallelChannels.PortHandle.LConn(2), 'autorouting', 'on') % Works

% links
add_line(subsys, UnitDelay.PortHandle.Outport,Pack.PortHandle.Inport, 'autorouting', 'on') % Works

% Storage:
% add_line(model, sourceLineHandles.Outport, destinationBlock, 'autorouting', 'on', 'arrows', 'on', 'srcPort', num2str(sourcePort), 'dstPort', num2str(destinationPort));
% add_line(model, CCS.LineHandle.LConn, Pack.LineHandle.LConn)
% 'autorouting', 'on', 'arrows', 'on', 'srcPort', num2str(sourcePort), 'dstPort', num2str(destinationPort));

%% Save results

save_system(model)
close_system(model)
battery_PackReplacement = 1; %Flag of function state   
fprintf ("Successfully replaced the BatteryPack\n")

%% __________________________________________________________________
% LossOutport_MA = add_block('nesl_utility/Probe', outport_MA1Name)
% add_line('Efficiency_Testbench/Pack1/ModuleAssembly1',)
%
% for i = 1:numel(connections)
%     sourceBlock = connections(i).SrcBlock;
%     sourcePort = connections(i).SrcPort;
%     fprintf('Connection %d: Source Block: %s, Source Port: %s\n', i, sourceBlock, sourcePort);
% end

%{
Cell_loss = out.bat_loss_out.Data % power_dissipation []
Bat_loss = Cell_loss * battery.ModInfo.num_cells_mod * battery.SysInfo.num_mods_sys;
Bat_out = abs(reshape(out.sim_out.battery_power_out.Data,[numel(out.sim_out.battery_power_out.Data), 1]));

Wirkungsgrad = Bat_out ./ (Bat_out + Bat_loss);
%}

%%_____________ dev help
% blockList = find_system(model, 'SearchDepth', 1, 'Type', 'Block');
% get_param(blockName,"MaskNames")
% get_param(blockName,"LineHandles")
% outputSrcBlock = get_param(lineHandles.LConn, 'SrcBlock')
% outputDstBlock = get_param(lineHandles.LConn, 'DstBlock')
% delete_line(lineHandles.LConn); % Works

end