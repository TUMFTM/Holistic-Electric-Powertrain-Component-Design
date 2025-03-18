function configs = append_configs_step_1(BatPara, BTMSPara, SysPara, SysSpec, configs)
    % Check if first run in this simulation. If first run, write at first array
    % position, otherwise append data to already existing data.

    persistent counter

    % Initialize the counter if it's the first run
    if isempty(counter)
        counter = 0;
    end

    % Check if the size of `configs` indicates a reset is needed
    if isstruct(configs) && numel(configs) == 1
        % Check if the `mod_ID` field exists and its value is NaN
        if isfield(configs, 'mod_ID') && isnan(configs.mod_ID)
           counter = 0;
        end   
    end 
    % Increment the counter to determine the next index
    counter = counter + 1;

    % Use the counter as the index for appending data
    configs(counter).mod_ID = counter;
    configs(counter).cell_ID = BatPara.name;
    configs(counter).BatPara = BatPara;
    configs(counter).BTMSPara = BTMSPara;
    configs(counter).SysPara = SysPara;
    configs(counter).SysSpec = SysSpec;
end