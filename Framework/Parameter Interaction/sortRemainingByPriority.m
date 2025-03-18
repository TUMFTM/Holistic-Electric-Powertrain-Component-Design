% Helper function to sort remaining variables by priority
function sortedVars = sortRemainingByPriority(remainingVars, priorityList)
    % Convert symbolic variables to strings
    remainingVarNames = arrayfun(@char, remainingVars, 'UniformOutput', false);
    
    % Sort the remaining variables based on the defined priority
    [~, sortedIndices] = sort(cellfun(@(x) find(strcmp(priorityList, x)), remainingVarNames));
    sortedVars = remainingVars(sortedIndices);
   
end