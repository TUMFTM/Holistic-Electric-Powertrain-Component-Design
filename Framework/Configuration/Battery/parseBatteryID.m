function config_battery = parseBatteryID(battery_GD, config_battery)
    % PARSEBATTERYID - Extracts cell type, chemistry, and capacity from battery ID.
    %
    % Inputs:
    %   - battery_GD: Struct containing the field 'cell_ID' (string) representing the battery name.
    %   - config_battery: Struct to store the parsed battery properties.
    %
    % Outputs:
    %   - config_battery: Updated struct with extracted cell type, chemistry, and capacity.

    % ---------------- Extract Cell Type ----------------
    if contains(battery_GD.cell_ID, "Pouch")
        config_battery.Cell_type = 2; % Pouch
    elseif contains(battery_GD.cell_ID, "Cyl")
        config_battery.Cell_type = 1; % Cylindrical
    elseif contains(battery_GD.cell_ID, "Pris")
        config_battery.Cell_type = 3; % Prismatic
    else
        error('Unknown cell type in battery name: %s', battery_GD.cell_ID);
    end

    % ---------------- Extract Cell Chemistry ----------------
    if contains(battery_GD.cell_ID, "NMC111")
        config_battery.Cell_chemistry = 1; % NMC111
    elseif contains(battery_GD.cell_ID, "NMC622")
        config_battery.Cell_chemistry = 2; % NMC622
    elseif contains(battery_GD.cell_ID, "NMC721")
        config_battery.Cell_chemistry = 3; % NMC721
    elseif contains(battery_GD.cell_ID, "NMC811")
        config_battery.Cell_chemistry = 4; % NMC811
    elseif contains(battery_GD.cell_ID, "LFP")
        config_battery.Cell_chemistry = 5; % LFP
    else
        error('Unknown cell chemistry in battery name: %s', battery_GD.cell_ID);
    end

    % ---------------- Extract Cell Capacity ----------------
    tokens = regexp(battery_GD.cell_ID, '_(\d+)Ah', 'tokens'); % Extract number before 'Ah'
    if ~isempty(tokens)
        config_battery.Cell_capacity = str2double(tokens{1}{1}); % Convert to numeric value
    else
        error('Cell capacity not found in battery name: %s', battery_GD.cell_ID);
    end
end