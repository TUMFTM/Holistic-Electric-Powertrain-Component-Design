function insertTableIntoWord(matlabTable, fileName)
    % Inserts a MATLAB table into a Word document.
    % Inputs:
    %   matlabTable - The MATLAB table to be inserted
    %   fileName - Name of the Word document to create (e.g., 'output.docx')

    % Start an ActiveX server for Word
    wordApp = actxserver('Word.Application');
    wordApp.Visible = true; % Set to true to see the Word window, false to run in the background

    % Create a new Word document
    doc = wordApp.Documents.Add;

    % Insert a title
    selection = wordApp.Selection;
    selection.TypeText('Exported Table from MATLAB');
    selection.TypeParagraph;

    % Get table dimensions
    [rows, cols] = size(matlabTable);

    % Add a table in Word
    wordTable = doc.Tables.Add(selection.Range, rows + 1, cols); % +1 for headers
    wordTable.Borders.Enable = true; % Add borders to the table

    % Insert headers into the table
    for col = 1:cols
        wordTable.Cell(1, col).Range.Text = matlabTable.Properties.VariableNames{col}; % Header text
        wordTable.Cell(1, col).Range.Bold = true; % Make headers bold
        wordTable.Cell(1, col).Range.ParagraphFormat.Alignment = 1; % Center align
    end

    % Insert data into the table
    for row = 1:rows
        for col = 1:cols
            % Convert numeric data to string, keep text as is
            if isnumeric(matlabTable{row, col})
                cellValue = num2str(matlabTable{row, col});
            else
                cellValue = matlabTable{row, col};
            end
            wordTable.Cell(row + 1, col).Range.Text = cellValue; % Data rows
        end
    end

    % Save the Word document
    fullPath = fullfile(pwd, fileName);
    doc.SaveAs(fullPath);
    disp(['Table inserted and saved to Word file: ', fullPath]);

    % Close Word
    doc.Close;
    wordApp.Quit;
    delete(wordApp);
end