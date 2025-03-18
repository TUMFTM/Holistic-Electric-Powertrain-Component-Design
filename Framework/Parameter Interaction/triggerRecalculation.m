    function globalParameters = triggerRecalculation(changedParam, fields, parameters, original_equations, priorityList, nonedit, constants, limits)

    global globalParameters;
    
    constants = evalin('base', 'constants');
    
    
    % Parameters where lower values are preferred
    % Define priority for each equation (lower number = higher priority)
    equationPriorities = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]; % Example priorities
    %equationPriorities = [6,7,8,9,10,11,12 ,1,2,3,4,5];
    % Sort equations by priority
    [~, sortedIndices] = sort(equationPriorities);
    original_equations = original_equations(sortedIndices); % Sort the equations list
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Substitute all constants in the equations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    equations = cell(size(original_equations));
    
    for i = 1:length(original_equations)
        equations{i} = subs(original_equations{i}, ...
            {'motor_specific_weight', 'gravitation', 'tire_coefficent', 'density', 'Cw', 'FrontalArea', ...
             'stages', 'const_mass', ...
             'motor_corner_speed'}, ...
            {constants.motor_specific_weight, constants.gravitation, constants.tire_coefficient, constants.density, ...
             constants.Cw, constants.FrontalArea, ...
             constants.stages, constants.const_mass, ...
             constants.motor_corner_speed});
    end
    %disp(constants.motor_specific_weight)
    %fprintf(constants.motor_specific_weight)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Substitute all constants in the equations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Store the current value of changedParam from globalParameters
    previousValue = globalParameters.(changedParam);
    fprintf('%s\n', changedParam);
    % Initialize solvedSet and recalculation counters for each variable
    solvedSet = {changedParam};
    recalcCount = struct();
    paramNames = fieldnames(parameters);
    for i = 1:length(paramNames)
        recalcCount.(paramNames{i}) = 0;  % Initialize counter to zero for each variable
    end
    
    % Set maximum recalculation attempts for each variable
    maxRecalcAttempts = 10;
    
    % Store recalculated values temporarily to validate them before updating
    tempParameters = globalParameters;
    %tempGlobalParameters = globalParameters;  % Store temporary global parameters

    % Loop through all fields and update tempParameters based on the field type
    % Update the parameter based on the changed field
    if isa(fields.(changedParam), 'matlab.ui.control.DropDown')
        parameters.(changedParam) = str2double(fields.(changedParam).Value);  % Convert dropdown value to double
        tempParameters.(changedParam) = parameters.(changedParam);  % Update global parameters
    else
        parameters.(changedParam) = fields.(changedParam).Value;  % Use numeric value directly
        tempParameters.(changedParam) = parameters.(changedParam);  % Update global parameters
    end
    
    %paramNames = fieldnames(fields);
    %for k = 1:length(paramNames)
    %    paramName = paramNames{k};
    %    
    %    % Update tempParameters with values from globalParameters if applicable
    %    if isfield(globalParameters, paramName) && ~isfield(fields, paramName)
    %        tempParameters.(paramName) = globalParameters.(paramName);
    %    end
    %end
    
    % Ensure all field values are reflected in tempParameters, converting dropdown values to numeric
    paramNames = fieldnames(fields);
    for k = 1:length(paramNames)
        paramName = paramNames{k};
        
        % Check if the field is a dropdown and convert to double if necessary
        if isa(fields.(paramName), 'matlab.ui.control.DropDown')
            tempParameters.(paramName) = str2double(fields.(paramName).Value);  % Convert dropdown value to double
        else
            tempParameters.(paramName) = fields.(paramName).Value;  % Use numeric value directly
        end
    end

    % Initialization part
    initialized = true;
    while initialized
        initialized = false;
        for i = 1:length(equations)
            
            eqn = equations{i};
            temp_eqn = eqn;
            varsInEqn = symvar(eqn);
            zeroVars = {};
            noneditiszero = false;

            for j = 1:length(varsInEqn)
                varName = char(varsInEqn(j));
                if isfield(tempParameters, varName)
                    currentVal = tempParameters.(varName);
                else
                        error('Unknown parameter');
                end

                if currentVal == 0
                    if ~isfield(nonedit, varName)
                        zeroVars{end+1} = varName;
                    else
                        noneditiszero = true;
                    end    
                end
            end

            for j = length(varsInEqn):-1:1
                varName = char(varsInEqn(j));
                if isfield(tempParameters, varName)
                    if tempParameters.(varName) ~= 0
                        eqn = subs(eqn, sym(varName), tempParameters.(varName));
                    end
                end
            end

            for j = 1:length(varsInEqn)
                varName = char(varsInEqn(j));
                if isfield(tempParameters, varName) && tempParameters.(varName) ~= 0
                    eqn = subs(eqn, sym(varName), tempParameters.(varName));
                end
            end

            if  length(zeroVars) == 1 && ~noneditiszero
                missingVar = zeroVars{1};
                disp(['The following variable will be initialized: ', varName]);
                %disp(['It is the only non-zero varible in the equation: ', char(temp_eqn)]);
                solution = solve(eqn, sym(missingVar));
                %solutions = real(solutions);
                %solution = double(solutions(1));

                if isempty(solution)
                    fprintf('No solution found for %s\n', missingVar);
                else
                    solvedValue = double(subs(solution, tempParameters));
                    if ~isempty(solvedValue)
                        tempParameters.(missingVar) = solvedValue;  % Temporarily store in tempGlobalParameters
                        initialized = true;
                    else
                        %fprintf('Error: Solved value for %s is empty\n', missingVar);
                    end
                end
            end
        end
    end

    % Recalculation part
    recalculationNeeded = true;
    while recalculationNeeded
        recalculationNeeded = false;
        for i = 1:length(equations)
            
            eqn = equations{i};
            tempo_eqn =eqn;
            varsInEqn = symvar(eqn);
            allSet = true;

            % Check if all variables in the equation have been set
            for j = 1:length(varsInEqn)
                varName = char(varsInEqn(j));
                if isfield(tempParameters, varName)
                    currentVal = tempParameters.(varName);
                else
                    allSet = false;
                    break;
                end
                if currentVal == 0
                    allSet = false;
                    break;
                end
            end

            % Substitute known values into the equation
            for j = length(varsInEqn):-1:1
                varName = char(varsInEqn(j));
                currentVal = tempParameters.(varName);
                if allSet && isfield(tempParameters, varName) && (~isfield(fields, varName) || isfield(nonedit, varName))
                    eqn = subs(eqn, sym(varName), tempParameters.(varName));
                    varsInEqn(j) = [];
                end
            end

            % Solve for remaining unknown variable if possible
            if allSet
                
                subParameters = tempParameters;
                fieldsToRemove = intersect(fieldnames(tempParameters), fieldnames(nonedit));
                
                subParameters = rmfield(subParameters, fieldsToRemove);
                lhsValue = double(subs(lhs(eqn), subParameters));
                rhsValue = double(subs(rhs(eqn), subParameters));

                tolerance = 1; % Set a suitable tolerance value
                if abs(lhsValue - rhsValue) > tolerance
                    %fprintf('A recalculation is necessary in the equation %s\n', eqn);
                    remainingVars = setdiff(varsInEqn, sym(changedParam));
                    remainingVars = setdiff(remainingVars, [sym(solvedSet), sym(changedParam)]);
                    remainingVars = sortRemainingByPriority(remainingVars, priorityList);
                    %remainingVars = sortRemainingByEquationPriority(remainingVars, eqn, equations, equationPriorities);
                    %if isempty(remainingVars)
                    %    fprintf('WARNING: no variable was found to be recalculated in this equation to keep it equal! The one with the least priority will just be recalculated')
                    %    remainingVars = setdiff(varsInEqn, sym(changedParam));
                    %    remainingVars = sortRemainingByPriority(remainingVars, priorityList);
                    %    solveForVar = char(remainingVars(end));
                    %end

                    if ~isempty(remainingVars)
                        solveForVar = char(remainingVars(end));
                        solution = solve(eqn, sym(solveForVar));
                        solvedValue = double(subs(solution, tempParameters));

                        if isreal(solvedValue) && ~isempty(solvedValue)
                            if isfield(tempParameters, solveForVar) && recalcCount.(solveForVar) > 0
                               
                                    fprintf('Skipping recalculation of %s because it was already recalculated.\n', solveForVar);
                                
                            else
                                % Assign the value for the first time or if recalcCount <= 0
                                tempParameters.(solveForVar) = solvedValue; % Temporarily store in tempGlobalParameters
                                solvedSet{end+1} = solveForVar;
                                recalculationNeeded = true;
                                recalcCount.(solveForVar) = recalcCount.(solveForVar) + 1;
                                break
                            end
                        else
                            fprintf('couldnt find a real solution')
                        end
                    else
                        disp('WARNING: no variable was found to be recalculated in this equation to keep it equal! Might be due to a circular dependency in the equations or making too many parameters uneditable')
                        disp(tempo_eqn)
                    end
                end
            end
        end
    end

    validationPassed = true;
    paramNames = fieldnames(tempParameters);  % Get all parameter names from tempParameters
    
    for k = 1:length(paramNames)
        paramName = paramNames{k};  % Extract the actual parameter name
        if isfield(limits, paramName)  % Check if the limits are defined for this parameter
            % Check if the parameter value is within the slider's limits
            if tempParameters.(paramName) < limits.(paramName)(1) || tempParameters.(paramName) > limits.(paramName)(2)
                fprintf('Validation failed: %s = %f is outside the allowed range [%f, %f]\n', ...
                    paramName, tempParameters.(paramName), limits.(paramName)(1), limits.(paramName)(2));
                validationPassed = false;
            end
        end
    end

    if validationPassed
        % Update parameters and globalParameters only if all values are valid
        parameters = tempParameters;
        globalParameters = tempParameters;  % Now update globalParameters

        % Update sliders and fields after validation
        for paramName = fieldnames(tempParameters)'
            paramName = char(paramName);
           
            %fprintf('Parameter: %s, Value: %f\n', paramName, tempParameters.(paramName));
         
            if isa(fields.(paramName), 'matlab.ui.control.Slider')
                fields.(paramName).Value = tempParameters.(paramName);  % Update slider value
            elseif isa(fields.(paramName), 'matlab.ui.control.NumericEditField')
                fields.(paramName).Value = tempParameters.(paramName);  % Update edit field value
            end
        end
    else
        fprintf('Recalculation aborted: One or more parameters are outside of the allowable limits.\n');
         % Reset the field value to the previous value
         
        if isa(fields.(changedParam), 'matlab.ui.control.DropDown')
            fields.(changedParam).Value = num2str(previousValue);  % Reset dropdown to the previous value
        else
            fields.(changedParam).Value = previousValue;  % Reset numeric field to the previous value
        end
    end

    % Loop through globalParameters and save each parameter to the base workspace
    paramNames = fieldnames(globalParameters);
    for k = 1:length(paramNames)
        paramName = paramNames{k};  % Get the parameter name
        paramValue = globalParameters.(paramName);  % Get the parameter value
        assignin('base', paramName, paramValue);  % Save to the base workspace
    end
end

function sortedVars = sortRemainingByEquationPriority(vars, currentEqn, allEquations, priorities)
    % Get the priority of the current equation
    eqnIndex = find(cellfun(@(x) isequal(x, currentEqn), allEquations));
    currentPriority = priorities(eqnIndex);

    % Keep only variables from equations with same or higher priority
    relevantVars = [];
    for i = 1:length(vars)
        for j = 1:length(allEquations)
            if ismember(vars(i), symvar(allEquations{j})) && priorities(j) <= currentPriority
                relevantVars = [relevantVars, vars(i)];
            end
        end
    end

    % Ensure unique values and sort them
    sortedVars = unique(relevantVars);
end
