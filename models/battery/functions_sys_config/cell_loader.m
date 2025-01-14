function [config_battery] = cell_loader(config_battery, chosen_capacity, bat_celltype, bat_cellchemistry)

%% cell_loader v0.2
% Translates cellnumber into loadable cell by adapting config_battery
% Part of Global Drive 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                           Versions                                    %
% v0.1 - first setup                                                    %
% v0.2 - special version to support battery Optimzation                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        known Problems                                 %
% -                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Normal NSGA2
if 0
    load("available_cells"); %#ok<LOAD>

    % adapt config_battery
    config_battery.cellname = available_cells.cellname(cell_idx_choosen);
    config_battery.manufacturer = available_cells.manufacturer(cell_idx_choosen);
    config_battery.manufacturer_parameterization = available_cells.manufacturer_parameterization(cell_idx_choosen);

    if available_cells.RC_class(cell_idx_choosen) ~=  "-"
        config_battery.simulation_charge_dynamics = convertStringsToChars(available_cells.RC_class(cell_idx_choosen));
        fprintf("battery: using %s charge dynamics based on cell availability \n",available_cells.RC_class(cell_idx_choosen))
    end
end

%% NSGA2_battery 
if 1
    % select data based on nsga input
    dataname = "available_cells";
    
    if bat_cellchemistry == 5
        bat_celltype = 3;
    end

    if bat_celltype == 1
        dataname = dataname + "_cyl";
        % Thin out for existing cells
        if ~ismember(cell_idx_choosen, [1, 5, 9])
            return
        end
    elseif bat_celltype == 2
        dataname = dataname + "_pouch";
    elseif bat_celltype == 3
        dataname = dataname + "_pris";
    else 
        return
    end
    
    if bat_cellchemistry == 1
        dataname = dataname + "_NMC111";
    elseif bat_cellchemistry == 2
        dataname = dataname + "_NMC622";
    elseif bat_cellchemistry == 3
        dataname = dataname + "_NMC721";
    elseif bat_cellchemistry == 4
        dataname = dataname + "_NMC811";
    elseif bat_cellchemistry == 5
        dataname = dataname + "_LFP";
    else
        return
    end
    
    l = load(dataname); %#ok<LOAD>

    % Extract capacities from cellname strings and find the closest match
    cellnames = l.(dataname).cellname;
    num_cells = numel(cellnames);
    capacities = zeros(1, num_cells);

    for i = 1:num_cells
        % Extract numeric value before "Ah" in the cellname
        capacity_match = regexp(cellnames{i}, '(\d+)(?=Ah)', 'match');
        if ~isempty(capacity_match)
            capacities(i) = str2double(capacity_match{1});
        else
            capacities(i) = NaN; % Handle cases where "Ah" is missing
        end
    end

    % Find the index of the closest capacity
    [~, closest_idx] = min(abs(capacities - chosen_capacity));

    % Set config_battery values based on the closest capacity
    config_battery.cellname = l.(dataname).cellname(closest_idx);
    config_battery.manufacturer = l.(dataname).manufacturer(closest_idx);
    config_battery.manufacturer_parameterization = l.(dataname).manufacturer_parameterization(closest_idx);
    
    if l.(dataname).RC_class(closest_idx) ~=  "-"
        config_battery.simulation_charge_dynamics = convertStringsToChars(l.(dataname).RC_class(closest_idx));
        fprintf("battery: using %s charge dynamics based on cell availability \n",l.(dataname).RC_class(closest_idx))
    end
end

%% Rest
% assignin('base','config_battery',config_battery)

clear available_cells
end