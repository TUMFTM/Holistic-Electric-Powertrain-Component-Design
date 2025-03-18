function function_handle_user_input(src, minVal, maxVal, paramName, fields, parameters, equations, priorityList, nonedit, outputFields, outputEquations, constants)
    % Check if the new value is within the specified limits
    if src.Value < minVal || src.Value > maxVal
        % Retrieve the last valid value from the workspace, defaulting to 0 if it does not exist
        if evalin('base', sprintf('exist(''%s'', ''var'')', paramName))
            lastValidValue = evalin('base', paramName);
        else
            lastValidValue = 0;
        end
    
        % Reset the field to the last valid value
        src.Value = lastValidValue;
    else
        % Trigger recalculation since the value is valid
        globalParameters = triggerRecalculation(paramName, fields, parameters, equations, priorityList, nonedit, outputEquations, constants);
        

        outputParameters = calculateOutputs(outputEquations, globalParameters, constants);
        
        updateOutputFields(outputParameters, fields);

        % After recalculation, update output fields with new values
        updateOutputs(outputFields);

    end
end