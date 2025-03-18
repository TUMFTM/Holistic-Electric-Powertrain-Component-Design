function [configs_4_sys_passed, configs_4_sys_failed] = calc_sys_dependent_variables(configs_3_sys_all, user)
    
    
    % calc_sys_dependent_variables: Evaluates system configurations against
    % predefined criteria and categorizes them into passed and failed configurations.
    %
    % Functionality:
    %   _ based on the selected number of parallel and serial
    %   connections various parameters on pack-level are calculated and
    %   saved for each configuration including: Pack dimesnions
    %   (without BTMS) and Pack Energy
    %   - Tests configurations for compliance with dimensional, energy, voltage, and current criteria.
    %   - Separates configurations into passed and failed categories.
    %   - Prepares the configurations for further processing.
    %
    % Inputs:
    %   - configs_3_sys_all: Cell array containing all system configurations.
    %   - user: Boolean flag to determine the inclusion of user-defined tests.
    %
    % Outputs:
    %   - configs_4_sys_passed: Cell array of configurations that passed all tests.
    %   - configs_4_sys_failed: Cell array of configurations that failed at least one test.

    % Preallocate arrays for passed and failed configurations
    configs_4_sys_passed = preallocate_configs_4_sys(); 
    configs_4_sys_failed = preallocate_configs_4_sys();

    % Iterate through all configurations
    for ii = 1:size(configs_3_sys_all, 2)
        config_bat = configs_3_sys_all(ii);

        % Test 1: Dimension constraints without BTMS (Battery Thermal Management System)
        [config_bat, passed_dim] = sys_test_dimensions_no_BTMS(config_bat);

        % Test 2: Energy content criteria
        [config_bat, passed_energy] = sys_test_energy(config_bat);

        % Test 3: Maximum voltage criteria
        [config_bat, passed_voltage] = sys_test_voltage(config_bat);

        % Test 4: Maximum current criteria
        [config_bat, passed_current] = sys_test_current(config_bat);

        % Combine individual test results into a single structure
        passed_sys1 = join_passed_structs(passed_dim, passed_energy);
        passed_sys2 = join_passed_structs(passed_voltage, passed_current);
        passed_sys = join_passed_structs(passed_sys1, passed_sys2);

        % Categorize configurations into passed or failed based on test results
        if check_for_failed_tests(passed_sys) && ~user
            % Configuration failed one or more tests
            configs_4_sys_failed = append_configs(configs_4_sys_failed, config_bat, passed_sys, 'fail');
        else
            % Configuration passed all tests
            configs_4_sys_passed = append_configs(configs_4_sys_passed, config_bat, passed_sys, 'pass');
        end
    end
end