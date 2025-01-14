function [] = displayStruct(myStruct,text1, text2, text3, text4, ignoreList)
% Function to output structs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = fieldnames(myStruct);
    for i = 1:numel(fields)
        field = fields{i};
        value = myStruct.(field);
        if ismember(field, ignoreList) ~= 1 
            if isnumeric(value) == 1
                disp([text1, field, text2, num2str(value), text3]);
            elseif istable(value)
                variable_names = value.Properties.VariableNames;
                variable_values= value.Variables;
                for j = 1: numel(variable_names)
                disp([text1, field, text2, char(variable_names(j)), text3, char(num2str(variable_values(j))), text4]);
                end
            else
                disp([text1, field, text2, char(value),text3]);
            end
        end
    end
end