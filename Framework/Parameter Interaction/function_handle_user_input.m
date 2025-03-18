function function_handle_user_input(src, minVal, maxVal, paramName, fields, parameters, equations, priorityList, nonedit, outputEquations, constants, limits)
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
        init_recalculate_update(paramName, fields, parameters, equations, priorityList, nonedit, outputEquations, constants, limits);
        % Trigger recalculation since the value is valid
        

    end
end