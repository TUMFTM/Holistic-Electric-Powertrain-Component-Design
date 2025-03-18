function cellTable = buildCellTable(cells)
    % Helper function to build a table from the cell data
    
    % Initialize arrays to store extracted data
    names = {};
    capacities = [];
    voltages = [];
    chemistries = {};
    dimensions = {};
    energies = [];
    currents = [];
    volumes = [];
    masses = [];
    specificEnergies = [];

    % Loop through each cell struct
    for i = 1:numel(cells)
        try
            % Extract data from the BatPara structure
            name = cells(i).Name;
            batPara = cells(i).BatPara;

            % Extract relevant parameters
            capacity = batPara.electrical.C_A; % Capacity in Ah
            voltage = batPara.electrical.U_nom; % Nominal voltage in V
            current = batPara.electrical.I_max; % Maximum current in A
            energy = batPara.electrical.E; % Energy in Wh
            chemistry = batPara.cellchemistry;
            dim_x = batPara.physical.dim_x*100; % Dimension x in m
            dim_y = batPara.physical.dim_y*100; % Dimension y in m
            dim_z = batPara.physical.dim_z*100; % Dimension z in m
            mass = batPara.physical.m; % Mass in kg

            % Determine the volume based on the cell type
            if contains(lower(name), 'cyl') % Cylindrical cells
                radius = dim_x / 2; % Assuming dim_x is the diameter
                volume = pi * (radius^2) * dim_y; % Volume of a cylinder
            else % Prismatic or pouch cells
                volume = dim_x * dim_y * dim_z; % Volume of a rectangular prism
            end

            % Calculate specific energy
            specificEnergy = energy / mass; % Specific energy in Wh/kg

            % Append data to arrays
            names{end + 1} = name; %#ok<AGROW>
            capacities(end + 1) = capacity; %#ok<AGROW>
            voltages(end + 1) = voltage; %#ok<AGROW>
            chemistries{end + 1} = chemistry; %#ok<AGROW>
            dimensions{end + 1} = sprintf('%.2f x %.2f x %.2f', dim_x, dim_y, dim_z); %#ok<AGROW>
            energies(end + 1) = energy; %#ok<AGROW>
            currents(end + 1) = current; %#ok<AGROW>
            volumes(end + 1) = volume; %#ok<AGROW>
            masses(end + 1) = mass; %#ok<AGROW>
            specificEnergies(end + 1) = specificEnergy; %#ok<AGROW>
        catch
            warning('Could not process cell: %s', cells(i).Name);
        end
    end

    % Build the table
    cellTable = table(names', capacities', voltages', chemistries', dimensions', energies', currents', volumes', masses', specificEnergies', ...
        'VariableNames', {'Name', 'Capacity (Ah)', 'Voltage (V)', 'Chemistry', 'Dimensions (mm)', 'Energy (Wh)', 'Current (A)', 'Volume (m^3)', 'Mass (kg)', 'Specific Energy (Wh/kg)'});
end