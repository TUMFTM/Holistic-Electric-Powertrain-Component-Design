function generateInteractivePlot(dropdownX, dropdownY, cellTable, cells)
    % Define plot options
    plotOptions = {'Capacity (Ah)', 'Voltage (V)', 'Energy (Wh)', 'Current (A)', 'Volume (m^3)', 'Mass (kg)', 'Specific Energy (Wh/kg)'};

    % Get selected options for X and Y axes
    selectedX = plotOptions{get(dropdownX, 'Value')};
    selectedY = plotOptions{get(dropdownY, 'Value')};

    try
        % Extract data for plotting
        xData = cellTable.(selectedX);
        yData = cellTable.(selectedY);

        % Filter out invalid data
        validIndices = ~isnan(xData) & ~isnan(yData) & ~isinf(xData) & ~isinf(yData);
        xData = xData(validIndices);
        yData = yData(validIndices);
        filteredCells = cells(validIndices); % Keep only valid cells

        % Create the plot
        figureHandle = figure('Name', 'Battery Cell Interactive Plot', ...
                              'NumberTitle', 'off', ...
                              'CloseRequestFcn', @(src, ~) delete(src), ...
                              'Position', [100, 100, 800, 600]);
        movegui(figureHandle, 'center');

        scatterHandle = scatter(xData, yData, 'filled');
        xlabel(selectedX, 'Interpreter', 'none');
        ylabel(selectedY, 'Interpreter', 'none');
        title('Battery Cell Data Plot', 'Interpreter', 'none');
        grid on;

        % Adjust axis limits for better visibility
        xlim([min(xData) - 5, max(xData) + 5]);
        ylim([min(yData) - 5, max(yData) + 5]);

        % Enable hover functionality
        dcm = datacursormode(figureHandle);
        dcm.Enable = 'on'; % Turn on data cursor mode
        dcm.DisplayStyle = 'datatip'; % Set hover style to 'datatip'

        % Inline logic to display only cell name, replacing underscores with spaces
        set(dcm, 'UpdateFcn', @(~, event) {
            ['Cell Name: ', strrep(filteredCells(event.DataIndex).Name, '_', ' ')]
        });

    catch ME
        % Handle errors and display useful information
        disp(['Error generating plot: ', ME.message]);
        if isfield(ME, 'stack')
            disp(ME.stack(1));
        end
    end
end