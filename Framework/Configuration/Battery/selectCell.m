function [battery_GD, BatPara] = selectCell(cells, battery_GD)
    % SELECTCELL - Opens a UI for the user to select a battery cell from a list.
    %
    % Inputs:
    %   cells      - Structure array containing available battery cells.
    %   battery_GD - Battery configuration structure to be updated.
    %
    % Outputs:
    %   battery_GD - Updated battery configuration with selected cell parameters.
    %   BatPara    - Battery parameters of the selected cell.
    BatPara = battery_GD.BatPara;
    
    % Get screen size for centering the UI figure
    screenSize = get(0, 'ScreenSize'); % Format: [x, y, width, height]
    screenWidth = screenSize(3);
    screenHeight = screenSize(4);

    % Define UI figure size
    figWidth = 400;
    figHeight = 200;

    % Calculate position to center the UI figure on the screen
    figX = (screenWidth - figWidth) / 2;
    figY = (screenHeight - figHeight) / 2;

    % Create the UI figure for cell selection
    fig = uifigure('Name', 'Select a Cell', ...
                   'Position', [figX, figY, figWidth, figHeight]);

    % Extract all available cell names from the structure
    cellNames = {cells.Name}; % Convert 'Name' field into a cell array

    % Create a dropdown menu populated with available cell names
    dropdown = uidropdown(fig, ...
        'Items', cellNames, ...                % Populate dropdown with cell names
        'Position', [100, 100, 200, 22], ...   % Set dropdown position and size
        'Value', cellNames{1});                % Set default selection to first item

    % Create a "Confirm Selection" button
    confirmButton = uibutton(fig, ...
        'Text', 'Confirm Selection', ...
        'Position', [150, 50, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) closeFigure()); % Assign callback function

    % Pause execution until user makes a selection and confirms
    
    disp('waiting')
    uiwait(fig);

    % Nested function to handle user confirmation and close the figure
    function closeFigure()
        selectedName = dropdown.Value; % Retrieve selected cell name

        % Find the corresponding cell structure in 'cells'
        selectedCell = cells(strcmp({cells.Name}, selectedName));

        % Update battery_GD with selected cell's properties
        battery_GD.BatPara = selectedCell.BatPara; % Assign selected cell parameters
        battery_GD.cell_ID = selectedCell.Name;    % Store selected cell name as ID
        BatPara = battery_GD.BatPara;              % Return selected cell parameters
        
        % Resume execution and close the UI
        uiresume(fig);
        pause(0.1)
        delete(fig);
    end
end