function battery_GD_All = process_selected_cells(selectedCells, dataTable, battery_GD, config_battery, config_vehicle)
    % Iterate through each row in selectedCells
    for i = 1:length(selectedCells)
   
        % Extract the selected cell name
        selectedCell = selectedCells(i);
        
        % Assign BatPara structure and cell name to battery_GD
        battery_GD.BatPara = selectedCell.BatPara;
        battery_GD.cell_ID = selectedCell.selectedCellName;
        
        config_battery = parseBatteryID(battery_GD, config_battery);

        [battery_GD, SysSpec, input_configs] = helper_initialization(battery_GD.BatPara, config_battery, battery_GD, config_vehicle);
 
        % Find the correct row in dataTable by comparing vehicle names
        rowIndex = find(strcmp(dataTable{:, 2}, selectedCell.vehicleName), 1); % Find first match
        
        if isempty(rowIndex)
            warning('No matching row found for vehicle: %s', selectedCell.vehicleName);
            continue; % Skip this iteration if no match is found
        end

        % Get the values from the respective columns in the matched row
        cellsInSeries = dataTable{rowIndex, 'CellsInSeriesInModule_Overview_'};
        cellsInParallel = dataTable{rowIndex, 'CellsInParallelInModule_Overview_'};
   
        if isnan(cellsInSeries) || isnan(cellsInParallel)
            continue;
        end   

        % Assign the values to battery_GD.SysPara
        battery_GD.SysPara.s_mod = cellsInSeries;
        battery_GD.SysPara.p_mod = cellsInParallel;
        
        battery_GDmodule = buildModuleWrapper(battery_GD, input_configs);
        battery_GDmodule = battery_GDmodule(1);

        battery_GD = battery_GDmodule; 

        % Set p_sys and s_sys from the respective columns
        battery_GD.SysInfo.num_parallel_mods_sys = dataTable{rowIndex, 'CellsInParallelInBatteryPack_Overview_'} / cellsInParallel;
        battery_GD.SysInfo.num_serial_mods_sys = dataTable{rowIndex, 'CellsInSeriesInBatteryPack_Overview_'} / cellsInSeries;
      
        disp(battery_GD.SysInfo.num_serial_mods_sys)
        battery_GD = buildPackWrapper(battery_GD);
        battery_GD = battery_GD(1); 
        % Display or log for debugging purposes
        if mod(battery_GD.SysInfo.num_serial_mods_sys, 1) ~= 0
            battery_GD.num_serial_rounded = true; % Number is not whole (e.g., 7.6)
         
        else
            battery_GD.num_serial_rounded = false; % Number is whole (e.g., 4, 23)
        end
        
        % Define the target directory
        % Define the existing target directory
        % Define the existing target directory
        currentDir = pwd; % Get the current MATLAB directory
        targetDir = fullfile(currentDir, 'AutoBattery', 'Vehicle Battery configuration');
        
        % Ensure the directory exists
        if ~exist(targetDir, 'dir')
            error('Target directory "Vehicle Battery configuration" does not exist.');
        end
        
        % Sanitize vehicle name to create a valid file name
        sanitizedVehicleName = regexprep(selectedCell.vehicleName, '[^\w]', '_'); % Replace non-alphanumeric characters with underscores
        
        % Define the file name
        fileName = sprintf('battery_GD_%s.mat', sanitizedVehicleName);
        
        % Full file path
        filePath = fullfile(targetDir, fileName);
        
        % Save the battery_GD variable
        save(filePath, 'battery_GD');
        disp(['Saved battery_GD for ', selectedCell.vehicleName, ' at ', filePath]);
        battery_GD_All(i).Vehicle = selectedCell.vehicleName;
        battery_GD_All(i).battery = battery_GD;
    end
end