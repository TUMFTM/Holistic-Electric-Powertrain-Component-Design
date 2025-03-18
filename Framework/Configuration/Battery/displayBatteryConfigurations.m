function  [battery_GD, f] = displayBatteryConfigurations(battery_GD)
    % displayBatteryConfigurations: Creates a UI for displaying battery configurations
    % and allows selection of a configuration.

    % Close any existing figure window
    close all;

    % Create the UI figure
    f = figure('Name', 'Battery Configurations', ...
               'NumberTitle', 'off', ...
               'Position', [200, 200, 1100, 500], ... 
               'MenuBar', 'none', 'ToolBar', 'none', ...
               'CloseRequestFcn', @(src, ~) uiresume(src)); % Resume execution when closed

    % Determine the number of configurations in `battery_GD`
    numConfigs = size(battery_GD, 2);

    % Prepare data for the table
    tableData = cell(numConfigs, 11); % Adjusted to include all columns
    for i = 1:numConfigs
        config = battery_GD(i);
        try
           % Populate table data
            tableData{i, 1} = i; % Index
            tableData{i, 2} = config.SysInfo.num_parallel_mods_sys; % Original p_sys
            tableData{i, 3} = config.SysInfo.num_serial_mods_sys; % Original s_sys
            tableData{i, 4} = config.SysPara.s_sys; % Total Serial (s_sys * s_mod)
            tableData{i, 5} = config.SysPara.p_sys; % Total Parallel (p_sys * p_mod)
            
            % New Columns for Module Parallel and Serial Connections
            tableData{i, 6} = config.SysPara.p_mod; % Parallel per Module
            tableData{i, 7} = config.SysPara.s_mod; % Serial per Module
            
            % Shifting previous columns to new positions
            tableData{i, 8} = sprintf('%.2f x %.2f x %.2f', ...
                                      config.SysInfo.dim_x_sys, ...
                                      config.SysInfo.dim_y_sys, ...
                                      config.SysInfo.dim_z_sys); % Dimensions
            tableData{i, 9} = config.SysInfo.mass_sys; % Mass
            tableData{i, 10} = config.SysInfo.C_sys; % Capacity
            tableData{i, 11} = config.SysInfo.E_sys; % Energy
            tableData{i, 12} = config.SysInfo.U_nom_sys; % Nominal Voltage
            tableData{i, 13} = config.SysInfo.I_max_sys; % Max Current
        catch
            % Handle missing or invalid data
            tableData{i, 1} = i; % Index
            tableData(i, 2:end) = {'Invalid or missing data'};
        end
    end

    % Define table column names
    columnNames = {'Index', 'Original p_sys', 'Original s_sys', ...
               'Total Serial (s_sys * s_mod)', 'Total Parallel (p_sys * p_mod)', ...
               'Parallel per Module (p_mod)', 'Serial per Module (s_mod)', ... % New Columns
               'Dimensions (m)', 'Mass (kg)', 'Capacity (Ah)', ...
               'Energy (Wh)', 'Nominal Voltage (V)', 'Max Current (A)'};

    % Create the table
    uitable(f, ...
            'Data', tableData, ...
            'ColumnName', columnNames, ...
            'Position', [50, 200, 1000, 200], ...
            'ColumnEditable', false, ...
            'FontSize', 10);

    % Dropdown menu for selecting a configuration index
    uicontrol(f, 'Style', 'text', ...
                 'String', 'Select Configuration Index:', ...
                 'Position', [50, 140, 200, 20], ...
                 'HorizontalAlignment', 'left', ...
                 'FontSize', 10);

    dropdown = uicontrol(f, 'Style', 'popupmenu', ...
                            'String', arrayfun(@num2str, 1:numConfigs, 'UniformOutput', false), ...
                            'Position', [250, 130, 150, 30], ...
                            'FontSize', 10);

    % Button to confirm selection
    uicontrol(f, 'Style', 'pushbutton', ...
                 'String', 'Select Configuration', ...
                 'Position', [450, 50, 150, 50], ...
                 'FontSize', 12, ...
                 'FontWeight', 'bold', ...
                 'Callback', @(~, ~) confirmSelection(dropdown.Value));

    % Button to select the optimal configuration
    uicontrol(f, 'Style', 'pushbutton', ...
             'String', 'Select Optimal Configuration', ...
             'Position', [650, 50, 150, 50], ...
             'FontSize', 12, ...
             'FontWeight', 'bold', ...
             'Callback', @selectOptimalConfiguration);

    % Block execution until GUI is closed
    uiwait(f);

    % Nested function to confirm user selection
    function confirmSelection(selectedIndex)
        battery_GD = battery_GD(:, selectedIndex); % Get the selected configuration
        uiresume(f); % Resume execution
        delete(f); % Close the figure
    end

    % Nested function to select the optimal configuration
    function selectOptimalConfiguration(~, ~)
        % Logic for selecting optimal configuration
        battery_GD = configs_elector(battery_GD);
        battery_GD = sustainability_estimator(battery_GD);
        uiresume(f); % Resume execution
        delete(f); % Close the figure
    end
end