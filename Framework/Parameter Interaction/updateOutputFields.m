function updateOutputFields(parameters, fields)
    % Update fields with the calculated output parameters
    paramNames = fieldnames(parameters);
    for i = 1:length(paramNames)
        paramName = paramNames{i};
        if isfield(fields, paramName)
            fields.(paramName).Value = parameters.(paramName);
        end
    end
end