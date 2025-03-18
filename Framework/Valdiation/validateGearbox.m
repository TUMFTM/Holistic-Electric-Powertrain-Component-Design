function simOut = validateGearbox(motormap, config_motor, simOut, battery_GD, config_gearbox)
    % VALIDATE_GEARBOX - Validates gearbox, motor, and battery operation.
    %
    % Inputs:
    %   - simOut: Structure containing simulation results, including required torques.
    %   - motormap: Structure containing motor's available speed and torque grids.
    %   - config_motor: Struct containing motor configuration details.
    %   - battery_GD: Struct containing battery parameters.
    %   - config_gearbox: Struct containing gearbox configuration details.
    %
    % Outputs:
    %   - Displays validation results for the gearbox, motor, and battery.

    % Calculate required gearbox torque with efficiency adjustment
    g_effi_safe = max(simOut.sim_out.g_effi.Data, 1);
    gearbox_torque_required = simOut.sim_out.WheelTorque.Data ./ (config_gearbox.iges * g_effi_safe);
    time = simOut.sim_out.WheelTorque.Time;  % Time vector

    % Extract allowed maximum torque from gearbox configuration
    allowed_gearbox_torque = config_gearbox.M_max;

    % Identify torque violations
    torque_violations = gearbox_torque_required > allowed_gearbox_torque;

    % Compute key torque statistics
    max_required_torque = max(gearbox_torque_required);
    spare_torque = allowed_gearbox_torque - max_required_torque;

    % Display validation results with structured output
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    disp('⚙️  **Gearbox Validation Results**');
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');

    if any(torque_violations)
        fprintf('⚠️  **WARNING**: Gearbox torque exceeded the allowed limit of %.2f Nm.\n', allowed_gearbox_torque);
        fprintf('🔹 Total violations: %d (%.2f%% of the cycle)\n', sum(torque_violations), ...
                sum(torque_violations) / length(torque_violations) * 100);
        fprintf('🔹 Maximum excess torque: %.2f Nm\n', max(gearbox_torque_required(torque_violations) - allowed_gearbox_torque));
    else
        disp('✅  Gearbox torque limits were NOT exceeded during the cycle.');
        fprintf('🔹 Spare torque margin: %.2f Nm\n', spare_torque);
        fprintf('🔹 Maximum required torque: %.2f Nm\n', max_required_torque);
    end

    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');

    % Define TUM color scheme
    tum_blue = [0, 101, 189] / 255;
    tum_dark_blue = [0, 79, 158] / 255;
    tum_light_blue = [100, 160, 200] / 255;
    tum_grey = [162, 162, 162] / 255;

    % Plot gearbox torque validation results
    figure('Name', 'Gearbox Torque Validation', 'Position', [100, 100, 1400, 700]);
    plot(time, gearbox_torque_required, '-', 'Color', tum_dark_blue, 'LineWidth', 2.5, 'DisplayName', 'Required Torque (Nm)');
    hold on;
    yline(allowed_gearbox_torque, '--', 'Color', tum_light_blue, 'LineWidth', 2.5, 'DisplayName', 'Allowed Torque (Nm)');

    % Highlight violations in the plot
    if any(torque_violations)
        scatter(time(torque_violations), gearbox_torque_required(torque_violations), 80, 'x', ...
                'MarkerEdgeColor', tum_light_blue, 'DisplayName', 'Violations', 'LineWidth', 2);
    end

    % Improve labels, title, and legend
    xlabel('Time (s)', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Torque (Nm)', 'FontSize', 20, 'FontWeight', 'bold');
    title('Gearbox Torque Validation', 'FontSize', 22, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 16);
    grid on;

    % Ensure axis font sizes update correctly
    set(gca, 'FontSize', 24, 'FontWeight', 'bold', 'GridColor', tum_grey);

    % If gearbox passes validation, suggest motor and battery revalidation
    if ~any(torque_violations)
        fprintf('\n🔄  Gearbox is sufficient. Now, the motor and battery should be rechecked.\n');
        % Uncomment if you want automatic validation of motor & battery:
        % validateMotor(motormap, config_motor, simOut);
        % validateBattery(simOut, battery_GD);
    end
end