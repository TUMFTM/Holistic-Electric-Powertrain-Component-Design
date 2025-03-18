function [configs_2_mod_passed, configs_2_mod_failed] = calc_mod_dependent_variables(configs_1_mod_all, user)
    % calc_mod_dependent_variables
    % Evaluates battery module configurations based on various criteria like maximum dimension limits
    % (excluding BTMS) and segregates them into passed and failed categories. The function also calculates
    % key module parameters such as capacity, nominal voltage, energy, maximum voltage,
    % minimum voltage, and maximum current.
    %
    % Inputs:
    %   configs_1_mod_all: Cell array containing initial module configurations.
    %   user: Boolean flag indicating whether user-defined constraints should be considered.
    %
    % Outputs:
    %   configs_2_mod_passed: Cell array containing configurations that passed all tests.
    %   configs_2_mod_failed: Cell array containing configurations that failed any test.

    % Preallocate cell arrays for passed and failed configurations
    configs_2_mod_passed = preallocate_configs_2_mod(); 
    configs_2_mod_failed = preallocate_configs_2_mod();

    % Step 1: Iterate through all initial configurations
    for ii = 1:size(configs_1_mod_all, 2)
        % Extract current configuration
        config_bat = configs_1_mod_all(ii);

        % Step 2: Test against maximum dimension limits (excluding BTMS)
        [config_bat, passed_dim] = mod_test_dimensions_no_BTMS(config_bat);

        % Step 3: Test against energy content criteria
        [config_bat, passed_energy] = mod_test_energy(config_bat);

        % Step 4: Combine results of individual tests into a unified structure
        passed_mod = join_passed_structs(passed_dim, passed_energy);

        % Step 5: Initialize thermal parameters (if required)
        %config_bat.thermal = struct();

        % Step 6: Classify configurations based on test results
        if check_for_failed_tests(passed_mod) && ~user
            % If any test failed and user flag is not set, classify as failed
            configs_2_mod_failed = append_configs(configs_2_mod_failed, config_bat, passed_mod, 'fail');
        else
            % Otherwise, classify as passed
            configs_2_mod_passed = append_configs(configs_2_mod_passed, config_bat, passed_mod, 'pass');
        end
    end
end