function availableCells = loadAllCells()
    % Load all available "dummy" cells and their BatPara structures from the given directory.
    %
    % Args:
    %   cellDir (string): The path to the directory containing the cell data.
    %
    % Returns:
    %   availableCells (struct array): A struct array with fields:
    %     - Name: The name of the cell
    %     - BatPara: The loaded BatPara structure for the cell

    % Get the list of cell names from the directory
    cellDir = "models/battery/input_and_parameters/01_cell_data/";
    cellFolders = dir(cellDir); % List contents of the directory
    cellNames = {cellFolders([cellFolders.isdir] & ~startsWith({cellFolders.name}, '.')).name}; % Get folder names, exclude '.' and '..'

    % Filter cell names to only include those starting with "dummy"
    dummyCellNames = cellNames(startsWith(cellNames, 'dummy'));

    % Initialize the output struct array
    availableCells = struct('Name', {}, 'BatPara', {});

    % Iterate over each "dummy" cell folder
    for i = 1:numel(dummyCellNames)
        cellName = dummyCellNames{i};

        % Construct the .mat file path
        matFilePath = fullfile(cellDir, cellName, [cellName, '.mat']);

        % Check if the .mat file exists
        if isfile(matFilePath)
            % Load the BatPara structure from the file
            loadedData = load(matFilePath);

            % Store the cell name and BatPara structure in the output struct
            availableCells(end + 1).Name = cellName; %#ok<AGROW>
            availableCells(end).BatPara = loadedData.BatPara; %#ok<AGROW>
        else
            % Display a warning if the .mat file is missing
            warning('File not found: %s', matFilePath);
        end
    end
end