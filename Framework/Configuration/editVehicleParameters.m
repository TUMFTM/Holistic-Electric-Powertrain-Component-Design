function config_vehicle = editVehicleParameters(config_vehicle)
    % EDITVEHICLEPARAMETERS - Allows users to edit vehicle parameters via a UI.
    %
    % Inputs:
    %   config_vehicle - Structure containing vehicle parameters.
    %
    % Outputs:
    %   config_vehicle - Updated structure with edited values.
    %
    % Description:
    %   This function creates a UI figure where users can edit specific vehicle 
    %   parameters grouped by categories. Changes are saved upon clicking the "Save" button.

    % Create the UI figure
    fig = uifigure('Name', 'Edit Vehicle Parameters', 'Position', [100, 100, 600, 400]);

    % Create a scrollable panel
    scrollPanel = uipanel(fig, ...
        'Scrollable', 'on', ...
        'Position', [10, 10, fig.Position(3) - 20, fig.Position(4) - 20]); % Panel fills most of the figure

    % Define groups and parameters manually
    groups = {
        'Performance', { 'maximum_speed', 'zero_to_corner_speed_in_s'};
        'Aerodynamics', {'cw_air', 'frontal_area'};
        'Dimensions', {'vehicle_dim_x_mm', 'vehicle_dim_y_mm', 'vehicle_dim_z_mm', 'tire_radius'};
        'Weight', {'mass_motor', 'mass_battery', 'mass_gearbox', 'mass_total', 'mass_module_factor', 'rot_mass_factor'}
    };

    % Calculate the required height for all groups
    rowHeight = 30; % Approximate height for each row in a table
    groupSpacing = 50; % Spacing for labels and tables
    totalHeight = 0;

    % Precompute totalHeight for all groups
    for i = 1:size(groups, 1)
        parameters = groups{i, 2};
        totalHeight = totalHeight + numel(parameters) * rowHeight + groupSpacing;
    end
    totalHeight = totalHeight + 50; % Add space for the Save button

    % Create a uipanel to contain the layout inside the scrollable panel
    contentPanel = uipanel(scrollPanel, ...
        'Position', [0, 0, scrollPanel.Position(3), totalHeight]); % Exceeds scrollPanel height

    % Create the main layout inside the content panel
    grid = uigridlayout(contentPanel, [size(groups, 1) * 2 + 1, 1], ...
        'RowSpacing', 10, 'Padding', [10 10 10 10]);

    % Store entries for Save functionality
    entries = struct();

    for i = 1:size(groups, 1)
        % Group label
        groupName = groups{i, 1};
        parameters = groups{i, 2};

        % Add label for the group
        groupLabel = uilabel(grid, ...
            'Text', groupName, ...
            'FontWeight', 'bold', ...
            'FontSize', 12, ...
            'HorizontalAlignment', 'center');
        groupLabel.Layout.Row = (i - 1) * 2 + 1;
        groupLabel.Layout.Column = 1;

        % Prepare data for the table
        paramNames = parameters';
        paramValues = cellfun(@(p) config_vehicle.(p), parameters, 'UniformOutput', false)';
        tableData = [paramNames, paramValues];

        % Calculate table height dynamically with constraints
        tableHeight = max(numel(parameters) * rowHeight, 60); % Set a minimum height of 60
        tableHeight = min(tableHeight, 200); % Optional: Set a maximum height of 200
        
        % Update row height in the grid layout
        grid.RowHeight{(i - 1) * 2 + 2} = tableHeight;

        % Create the table
        t = uitable(grid, ...
            'Data', tableData, ...
            'ColumnName', {'Parameter', 'Value'}, ...
            'ColumnEditable', [false, true], ... % Only "Value" column is editable
            'RowName', {}, ... % Remove row numbering
            'ColumnWidth', {'1x', '1x'}); % Automatically adjust width
        t.Layout.Row = (i - 1) * 2 + 2; % Place below the label
        t.Layout.Column = 1;

        % Save references to update entries
        for j = 1:numel(parameters)
            paramName = parameters{j};
            entries.(paramName) = t;
        end
    end

    % Add Save button
    saveButton = uibutton(grid, ...
        'Text', 'Save', ...
        'ButtonPushedFcn', @(~, ~) saveChanges(), ...
        'HorizontalAlignment', 'center');
    saveButton.Layout.Row = size(groups, 1) * 2 + 1;
    saveButton.Layout.Column = 1;

    % Wait for user to close the figure
    waitfor(fig);

    % Save function
    function saveChanges()
        % Update config_vehicle with edited values
        for i = 1:size(groups, 1)
            groupParams = groups{i, 2};
            table = entries.(groupParams{1}); % Get the group's table
            editedData = table.Data;

            for j = 1:numel(groupParams)
                config_vehicle.(groupParams{j}) = editedData{j, 2};
            end
        end
        close(fig);
    end
end