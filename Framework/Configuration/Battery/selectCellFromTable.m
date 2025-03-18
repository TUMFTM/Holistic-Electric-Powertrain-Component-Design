function selectCellFromTable(event)
    % Handle selection of a cell in the table view and save its parameters as BatPara.
    %
    % Args:
    %   event: Cell selection event data.
    %   cells (struct array): Struct array containing detailed cell data.

    if isempty(event.Indices)
        return;
    end

   % Get the name of the selected cell
    selectedName = event.Source.Data{event.Indices(1), 1};
    
    % Feedback to the user
    fprintf('The cell clicked on is: %s\n', selectedName);
end