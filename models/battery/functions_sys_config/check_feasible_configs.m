function check_feasible_configs(configs,configs_2)

num_of_configs = size(configs, 2);

if isnan(configs(1).mod_ID)
    if isstruct(configs_2)
        %% Error message if all configs failed
        h.fnames = fieldnames(configs_2(1).Tests_mod);
        
        err_msg = ['\n########################################################', ...
            '\nNo feasible concepts in %s! \n', '\nThese configurations failed:'];
        err_msg = sprintf(err_msg, inputname(1));
        
        for i = 1:numel(configs_2)
            % Adapt Config ID
            if isfield(configs_2, 'Tests_BTMS')
                h.IDtype = 'sys_ID';
                t.Testtype = 'Tests_BTMS';
            elseif isfield(configs_2, 'Tests_sys')
                h.IDtype = 'sys_ID';
                t.Testtype = 'Tests_sys';
            elseif isfield(configs_2, 'Tests_mod')
                    h.IDtype = 'mod_ID';
                    t.Testtype = 'Tests_mod';
            end
            h.fnames = fieldnames(configs_2(1).(t.Testtype));

            err_msg_config = sprintf('\n\nConfiguration %i', configs_2(i).(h.IDtype));
            err_msg_var = {};
            
            for j = 1:numel(h.fnames)
                % Field name and value
                field_name = h.fnames{j};
                field_value = configs_2(i).(t.Testtype).(field_name);
                err_msg_field = sprintf('\n %s: %i', field_name, field_value);
                err_msg_var = [err_msg_var, {err_msg_field}];
            end
            
            % Concatenate the config ID and field name/value strings
            err_msg_config = [err_msg_config, err_msg_var{:}];
            err_msg = [err_msg, err_msg_config];
        end
        
        err_msg = [err_msg, '\n########################################################'];

        % Format the error message string with the input name and concatenated error message
        err_msg = sprintf(err_msg, inputname(1));
        
        % Throw the error with the formatted error message string
        error(err_msg);

    else
        error('No feasible concepts in %s!', inputname(1))
    end
else
   fprintf('There are %i concepts in %s.\n', num_of_configs ,inputname(1))
end