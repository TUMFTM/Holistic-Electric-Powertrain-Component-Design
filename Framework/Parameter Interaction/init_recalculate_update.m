function init_recalculate_update(paramName, fields, parameters, equations, priorityList, nonedit, outputEquations, constants, limits)
    % Trigger recalculation of global parameters based on the given input parameters
    globalParameters = triggerRecalculation(paramName, fields, parameters, equations, priorityList, nonedit, constants, limits);
    
    % Compute output parameters using the provided equations and global parameters
    outputParameters = calculateOutputs(outputEquations, globalParameters, constants);
    
    % Update the relevant output fields with the newly calculated parameters
    updateOutputFields(outputParameters, fields);
end