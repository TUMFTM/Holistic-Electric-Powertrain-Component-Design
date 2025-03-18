function showPlotView(cellTable, cells, parentFigure)
    % Display the plot view in the GUI with extended plot options.
    %
    % Args:
    %   cellTable (table): Table containing battery cell data.
    %   cells (struct array): Struct array containing detailed cell data.
    %   parentFigure (figure): The parent figure.

    % Clear existing content in the figure
    clf(parentFigure);

    % Dropdowns for selecting X and Y axes
    uicontrol('Style', 'text', 'String', 'X-axis:', ...
              'Position', [50, 500, 50, 25], 'HorizontalAlignment', 'left');
    plotOptions = {'Capacity (Ah)', 'Voltage (V)', 'Energy (Wh)', 'Current (A)', 'Volume (m^3)', 'Mass (kg)', 'Specific Energy (Wh/kg)'};
    dropdownX = uicontrol('Style', 'popupmenu', 'String', plotOptions, ...
                          'Position', [100, 500, 150, 25]);

    uicontrol('Style', 'text', 'String', 'Y-axis:', ...
              'Position', [300, 500, 50, 25], 'HorizontalAlignment', 'left');
    dropdownY = uicontrol('Style', 'popupmenu', 'String', plotOptions, ...
                          'Position', [350, 500, 150, 25]);

    % Button to generate the plot
    uicontrol('Style', 'pushbutton', 'String', 'Generate Plot', ...
              'Position', [550, 500, 100, 30], ...
              'Callback', @(~, ~) generateInteractivePlot(dropdownX, dropdownY, cellTable, cells));

    % Add navigation buttons again
    uicontrol('Style', 'pushbutton', 'String', 'Table View', ...
              'Position', [50, 550, 100, 30], ...
              'Callback', @(~, ~) showTableView(cellTable, cells));

    uicontrol('Style', 'pushbutton', 'String', 'Plot View', ...
              'Position', [200, 550, 100, 30], ...
              'Callback', @(~, ~) showPlotView(cellTable, cells, parentFigure));
end