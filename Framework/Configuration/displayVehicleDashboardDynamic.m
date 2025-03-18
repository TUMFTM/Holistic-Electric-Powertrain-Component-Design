function displayVehicleDashboardDynamic(config_vehicle, config_motor, battery_GD, config_gearbox, config_Inverter, motormap_name)
    % Create or reuse a persistent figure for the dashboard
    persistent dashFig timerObj;
    if isempty(dashFig) || ~isvalid(dashFig)
        dashFig = uifigure('Name', 'Vehicle Dashboard', 'Position', [100, 100, 500, 700]);

        % Set up a timer to update the dashboard every second
        timerObj = timer('ExecutionMode', 'fixedRate', 'Period', 1, ...
                         'TimerFcn', @(~, ~) updateDashboard(dashFig, config_vehicle, config_motor, battery_GD, config_gearbox, config_Inverter, motormap_name));
        start(timerObj);
    end

    % Initialize the dashboard layout
    updateDashboard(dashFig, config_vehicle, config_motor, battery_GD, config_gearbox, config_Inverter, motormap_name);
end

function updateDashboard(dashFig, config_vehicle, config_motor, battery_GD, config_gearbox, config_Inverter, motormap_name)
    % Clear the figure to refresh displayed values
    delete(dashFig.Children);

    % Title
    uilabel(dashFig, 'Text', 'Vehicle Information Dashboard', ...
        'Position', [100, 650, 300, 30], 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % General Vehicle Information
    addSection(dashFig, 'General Vehicle Information', [50, 600, 400, 20]);
    addLabel(dashFig, sprintf('Total Mass: %.2f kg', getFieldOrDefault(config_vehicle, 'mass_total', 0)), [50, 580, 400, 20]);
    addLabel(dashFig, sprintf('Inverter Type: %s', getFieldOrDefault(config_vehicle, 'inverter_type', 'N/A')), [50, 560, 400, 20]);

    % Battery Information
    addSection(dashFig, 'Battery Information', [50, 520, 400, 20]);
    addLabel(dashFig, sprintf('Mass Battery: %.2f kg', getFieldOrDefault(config_vehicle, 'mass_battery', 0)), [50, 500, 400, 20]);
    addLabel(dashFig, sprintf('Battery C_Sys: %.2f', getFieldOrDefault(battery_GD.SysInfo, 'C_sys', 0)), [50, 480, 400, 20]);
    addLabel(dashFig, sprintf('Nominal Voltage: %.2f V', getFieldOrDefault(battery_GD.SysInfo, 'U_nom_sys', 0)), [50, 460, 400, 20]);
    addLabel(dashFig, sprintf('Cell ID: %s', getFieldOrDefault(battery_GD, 'cell_ID', 'N/A')), [50, 440, 400, 20]);
    addLabel(dashFig, sprintf('Battery Energy: %.2f kWh', getFieldOrDefault(battery_GD, 'E_sys', 0)), [50, 420, 400, 20]);

    % Motor Information
    addSection(dashFig, 'Motor Information', [50, 380, 400, 20]);
    addLabel(dashFig, sprintf('Mass Motor: %.2f kg', getFieldOrDefault(config_vehicle, 'mass_motor', 0)), [50, 360, 400, 20]);
    addLabel(dashFig, sprintf('Peak Torque: %.2f Nm', getFieldOrDefault(config_motor, 'peak_torque', 0)), [50, 340, 400, 20]);
    addLabel(dashFig, sprintf('Peak Power: %.2f kW', getFieldOrDefault(config_motor, 'peak_power', 0)), [50, 320, 400, 20]);
    addLabel(dashFig, sprintf('Motor Name: %s', getFieldOrDefault(motormap_name, '', 'N/A')), [50, 300, 400, 20]);

    % Gearbox Information
    addSection(dashFig, 'Gearbox Information', [50, 260, 400, 20]);
    addLabel(dashFig, sprintf('Mass Gearbox: %.2f kg', getFieldOrDefault(config_vehicle, 'mass_gearbox', 0)), [50, 240, 400, 20]);
    addLabel(dashFig, sprintf('Gear Ratio (iges): %.2f', getFieldOrDefault(config_gearbox, 'iges', 0)), [50, 220, 400, 20]);
    addLabel(dashFig, sprintf('Maximum Torque: %.2f Nm', getFieldOrDefault(config_gearbox, 'M_max', 0)), [50, 200, 400, 20]);

    % Inverter Information
    addSection(dashFig, 'Inverter Information', [50, 160, 400, 20]);
    addLabel(dashFig, sprintf('Switching Frequency: %.2f Hz', getFieldOrDefault(config_Inverter, 'fswitch', 0)), [50, 140, 400, 20]);
end

function addSection(parent, sectionTitle, position)
    % Add a bold section title
    uilabel(parent, 'Text', ['--- ' sectionTitle ' ---'], ...
        'Position', position, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
end

function addLabel(parent, labelText, position)
    % Add a single-line label
    uilabel(parent, 'Text', labelText, ...
        'Position', position, 'HorizontalAlignment', 'left');
end

function value = getFieldOrDefault(structure, fieldName, defaultValue)
    % Safely get a field value or use a default
    if isstruct(structure) && isfield(structure, fieldName)
        value = structure.(fieldName);
    else
        value = defaultValue;
    end
end