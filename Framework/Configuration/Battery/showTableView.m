function showTableView(cellTable, cells)
    % Display the table view in the GUI with sorting options.
    %
    % Args:
    %   cellTable (table): Table containing battery cell data.
    %   cells (struct array): Struct array containing detailed cell data.
    %   parentFigure (figure): The parent figure.

    % Clear existing content in the figure
    % Create the main figure
    parentFigure = figure('Name', 'Battery Cell Viewer', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);

    % Create the table
    t = uitable('Parent', parentFigure, ...
                'Data', table2cell(cellTable), ...
                'ColumnName', cellTable.Properties.VariableNames, ...
                'Position', [50, 50, 700, 400], ... % Adjusted table position and size
                'CellSelectionCallback', @(src, event) selectCellFromTable(event, cells));

    % Create navigation buttons
    uicontrol('Style', 'pushbutton', 'String', 'Table View', ...
              'Position', [50, 500, 100, 30], ...
              'Callback', @(~, ~) showTableView(cellTable, cells, parentFigure));

    uicontrol('Style', 'pushbutton', 'String', 'Plot View', ...
              'Position', [200, 500, 100, 30], ...
              'Callback', @(~, ~) showPlotView(cellTable, cells, parentFigure));

    % Create dropdown for sorting options
    uicontrol('Style', 'text', 'String', 'Sort by:', ...
              'Position', [50, 460, 50, 25], 'HorizontalAlignment', 'left'); % Adjusted position

    sortOptions = {'Name', 'Capacity (Ah)', 'Voltage (V)', 'Chemistry', 'Dimensions (mm)'};
    dropdown = uicontrol('Style', 'popupmenu', 'String', sortOptions, ...
                         'Position', [110, 460, 150, 25]); % Adjusted position

    % Create a button to apply sorting
    uicontrol('Style', 'pushbutton', 'String', 'Sort', ...
              'Position', [280, 460, 50, 25], ... % Adjusted position
              'Callback', @(~, ~) sortTable(dropdown, t, cellTable));
end