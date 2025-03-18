function selectCellOnClick(~, event, cells, xData, yData)
    % Handle click events to select a cell and save its full BatPara struct.
    %
    % Args:
    %   ~: The graphics object triggering the event (not used here).
    %   event: The click event object.
    %   cells: The struct array containing detailed cell data (including BatPara).
    %   xData, yData: The data arrays used in the plot.
    
    % Get the clicked point's coordinates
    clickedPoint = event.IntersectionPoint(1:2);

    % Find the closest point in the data
    distances = vecnorm([xData - clickedPoint(1), yData - clickedPoint(2)], 2, 2);
    [~, idx] = min(distances);

    % Retrieve the full BatPara struct of the corresponding cell
    selectedCell = cells(idx);

    selectedName = selectedCell.Name;
    
    fprintf('The cell clicked on is: %s\n', selectedName);
end