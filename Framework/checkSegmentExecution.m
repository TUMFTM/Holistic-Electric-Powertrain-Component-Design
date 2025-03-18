function status = checkSegmentExecution(currentSegment)
    % Define dependencies between segments
    persistent segmentDependencies segmentStatus;

    % Initialize dependencies only once
    if isempty(segmentDependencies)
        segmentDependencies.Initialization = {};  
        segmentDependencies.Configuration = {{'Initialization'}};  
        segmentDependencies.EditConfiguration = {{'Configuration'}};  
        segmentDependencies.MotorSegmentInfo = {{'Configuration'}};  
        segmentDependencies.Motor1Selection = {{'Initialization'}};  
        segmentDependencies.Motor1Topology = {{'Initialization'}};  
        segmentDependencies.NoMotor = {{'Motor1Selection'}, {'Motor1Topology'}, {'Configuration'}};  
        segmentDependencies.Motor2Selection = {{'Initialization'}};  
        segmentDependencies.Motor2Topology = {{'Motor1Selection', 'Motor1Topology'}, {'Motor2Selection'}, {'Configuration'}};  
        segmentDependencies.Motor1Output = {{'Motor1Selection'}, {'Motor1Topology'}};  
        segmentDependencies.BatterySegmentInfo = {{'Configuration'}, {'Motor2Topology', 'NoMotor'}, {'Initialization'}};  
        segmentDependencies.BatteryMotorInfo = {{'Motor2Topology', 'NoMotor'}, {'Initialization'}};  
        segmentDependencies.cellDisplay = {'Initialization'};  
        segmentDependencies.cellSelection = {{'Initialization'}, {'Configuration'}}; 
        segmentDependencies.Module = {'cellSelection'}; 
        segmentDependencies.ModuleOutput = {'Module'}; 
        segmentDependencies.Pack = {'Module'}; 
        segmentDependencies.PackOutput = {'Pack'}; 
        segmentDependencies.GearboxDesign = {'Configuration'}; 
        segmentDependencies.GearboxOutput = {'GearboxDesign'}; 
        segmentDependencies.Inverter = {{'Motor2Topology', 'NoMotor'}, {'Configuration'}};
        segmentDependencies.ValidationSettings = {{'Initialization'}}; 
        segmentDependencies.ValidateMotor = {{'ValidationSettings'}, {'Motor2Topology', 'NoMotor'}}; 
        segmentDependencies.ValidateBattery = {{'ValidationSettings'}, {'ValidateMotor','ValidateGearbox'}, {'Pack'}}; 
        segmentDependencies.ValidateGearbox = {{'ValidationSettings'},  {'Motor2Topology', 'NoMotor'}, {'GearboxDesign'}}; 
        segmentDependencies.AttributesSettings = {{'Initialization'}}; 
        segmentDependencies.Performance = {{'AttributesSettings'}, {'Motor2Topology', 'NoMotor'}, {'Pack'}, {'GearboxDesign'}, {'Inverter'}};
        segmentDependencies.Consumption = { {'Motor2Topology', 'NoMotor'}, {'Pack'}, {'GearboxDesign'}, {'Inverter'}};
        segmentDependencies.Simulation = {{'Motor2Topology', 'NoMotor'}, {'Pack'}, {'GearboxDesign'}, {'Inverter'}};
        segmentDependencies.loadVehicle = {'Initialization'}; 
        segmentDependencies.loadBattery = {'Initialization'}; 
    end

    % Initialize persistent segment execution status if empty
    if isempty(segmentStatus)
        segmentStatus = struct();  
    end

    % Ensure segment exists
    if ~isfield(segmentDependencies, currentSegment)
        error('Unknown segment: %s', currentSegment);
    end
    
    % Special Case: If `loadVehicle` is executed, mark its implied segments as executed
    if strcmp(currentSegment, 'loadVehicle')
        impliedSegments = {'Configuration', 'Motor2Topology', 'Pack', 'GearboxDesign', 'Inverter'};
        for i = 1:numel(impliedSegments)
            segmentStatus.(impliedSegments{i}) = true;
        end
    end

    % Special Case: If `loadBattery` is executed, mark its implied segments as executed
    if strcmp(currentSegment, 'loadBattery')
        impliedSegments = {'cellSelection', 'Module', 'Pack'};
        for i = 1:numel(impliedSegments)
            segmentStatus.(impliedSegments{i}) = true;
        end
    end

    % Get dependencies for the current segment (list of alternative sets)
    dependencyGroups = segmentDependencies.(currentSegment);

    % If no dependencies, mark segment as executed
    if isempty(dependencyGroups)
        segmentStatus.(currentSegment) = true;
        status = true;
        return;
    end

    % Check if at least **one condition from each group** is satisfied
    allGroupsSatisfied = true;
    missingGroups = {};

    for i = 1:numel(dependencyGroups)
        dependencies = dependencyGroups{i}; % One group of possible dependencies

        % Ensure dependencies are a cell array of character vectors
        dependencies = cellstr(dependencies);  

        % Check if any dependency in this group is fulfilled
        isGroupSatisfied = any(isfield(segmentStatus, dependencies) & structfun(@(x) x, segmentStatus, 'UniformOutput', true));

        if ~isGroupSatisfied
            allGroupsSatisfied = false;
            missingGroups{end+1} = strjoin(dependencies, ' OR '); % Store missing group
        end
    end

    % If any group was not satisfied, display a warning and block execution
    if ~allGroupsSatisfied
        warning(['Missing dependencies: [' strjoin(missingGroups, '] AND [') ']']);
        status = false;
        return;
    end

    % If all dependencies are met, allow execution
    segmentStatus.(currentSegment) = true;
    status = true;
end