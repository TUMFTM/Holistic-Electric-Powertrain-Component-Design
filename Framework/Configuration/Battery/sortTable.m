function sortTable(dropdown, t, cellTable)
    % Sort the table based on the selected parameter
    %
    % Args:
    %   dropdown: The dropdown menu for selecting sorting parameter.
    %   t: The uitable to update.
    %   cellTable: The original table.

    % Get the selected sorting parameter
    sortOptions = {'Name', 'Capacity (Ah)', 'Voltage (V)', 'Chemistry', 'Dimensions (mm)'};
    selectedOption = sortOptions{get(dropdown, 'Value')};

    % Sort the table
    sortedTable = sortrows(cellTable, selectedOption);

    % Update the uitable
    t.Data = table2cell(sortedTable);
end