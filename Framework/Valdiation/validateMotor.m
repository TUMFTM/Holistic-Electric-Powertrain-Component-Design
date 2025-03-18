function validateMotor(motormap, config_motor, simOut)
    % VALIDATE_MOTOR - Validate motor speed and torque requirements
    % against the motor's capability defined in the motor map.
    %
    % Outputs:
    %   - Displays validation results and visualizations.

    % Extract required speed and torque from simulation output
    required_speed = simOut.sim_out.Motor_rpm.Data;  
    required_torque = simOut.sim_out.Motor_Torque.Data;  
    time = simOut.sim_out.Motor_rpm.Time;  

    % Extract motor map data
    motor_speed_map = motormap.Speed;  
    motor_torque_map = motormap.Shaft_Torque;

    % Initialize tracking variables
    numPoints = length(required_speed);
    violations_motor1 = false(numPoints, 1);
    violations_motor2 = false(numPoints, 1);
    available_torque_motor1 = zeros(numPoints, 1);
    available_torque_motor2 = zeros(numPoints, 1);
    max_deviation_motor1 = 0;
    max_deviation_motor2 = 0;
    max_spare_torque_motor1 = 0;
    max_spare_torque_motor2 = 0;

    % Determine torque split for dual-motor configurations
    if config_motor.torque2 > 0
        torque_ratio_motor1 = config_motor.power1 / (config_motor.power1 + config_motor.power2);
        torque_ratio_motor2 = config_motor.power2 / (config_motor.power1 + config_motor.power2);
        required_torque_motor1 = required_torque * torque_ratio_motor1;
        required_torque_motor2 = required_torque * torque_ratio_motor2;
    else
        required_torque_motor1 = required_torque;
        required_torque_motor2 = zeros(size(required_torque));
    end

    % Validate torque-speed pairs for each time step
    for i = 1:numPoints
        % Find the closest speed in the motor map
        [~, speed_row_idx] = min(abs(motor_speed_map(:, 1) - required_speed(i)));
        available_torques_motor1 = motor_torque_map(speed_row_idx, :) * config_motor.motorsCount1;

        % Validate Motor 1
        if required_torque_motor1(i) >= 0
            available_torque_motor1(i) = max(available_torques_motor1);
            if required_torque_motor1(i) > available_torque_motor1(i)
                violations_motor1(i) = true;
                max_deviation_motor1 = max(max_deviation_motor1, abs(required_torque_motor1(i) - available_torque_motor1(i)));
            else
                max_spare_torque_motor1 = max(max_spare_torque_motor1, available_torque_motor1(i) - required_torque_motor1(i));
            end
        end

        % Validate Motor 2 if applicable
        if config_motor.torque2 > 0
            available_torques_motor2 = motor_torque_map(speed_row_idx, :) * config_motor.motorsCount2;

            if required_torque_motor2(i) >= 0
                available_torque_motor2(i) = max(available_torques_motor2);
                if required_torque_motor2(i) > available_torque_motor2(i)
                    violations_motor2(i) = true;
                    max_deviation_motor2 = max(max_deviation_motor2, abs(required_torque_motor2(i) - available_torque_motor2(i)));
                else
                    max_spare_torque_motor2 = max(max_spare_torque_motor2, available_torque_motor2(i) - required_torque_motor2(i));
                end
            end
        end
    end

    % Display validation results with structured output
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    disp('⚙️  **Motor Validation Results**');
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');

    % Display results for Motor 1
    disp('🔹 **Motor 1:**');
    if any(violations_motor1)
        fprintf('⚠️  Motor 1 has %d violations.\n', sum(violations_motor1));
        fprintf('🔹 Maximum Torque Deviation: %.2f Nm\n', max_deviation_motor1);
    else
        disp('✅  Motor 1 operates within safe limits.');
        fprintf('🔹 Maximum Spare Torque: %.2f Nm\n', max_spare_torque_motor1);
    end

    % Display results for Motor 2 (if applicable)
    if config_motor.torque2 > 0
        disp('🔹 **Motor 2:**');
        if any(violations_motor2)
            fprintf('⚠️  Motor 2 has %d violations.\n', sum(violations_motor2));
            fprintf('🔹 Maximum Torque Deviation: %.2f Nm\n', max_deviation_motor2);
        else
            disp('✅  Motor 2 operates within safe limits.');
            fprintf('🔹 Maximum Spare Torque: %.2f Nm\n', max_spare_torque_motor2);
        end
    end

    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');

    % Define TUM colors
    tum_blue = [0, 101, 189] / 255;
    tum_dark_blue = [0, 79, 158] / 255;
    tum_light_blue = [100, 160, 200] / 255;
    tum_grey = [162, 162, 162] / 255;

    % Plot results for Motor 1
    figure('Position', [100, 100, 1400, 700]);
    plot(time, required_torque_motor1, '-', 'Color', tum_dark_blue, 'LineWidth', 2.5, 'DisplayName', 'Required Torque Motor 1');
    hold on;
    plot(time, available_torque_motor1, '--', 'Color', tum_light_blue, 'LineWidth', 2.5, 'DisplayName', 'Available Torque Motor 1');
    scatter(time(violations_motor1), required_torque_motor1(violations_motor1), 80, 'x', 'MarkerEdgeColor', tum_light_blue, 'DisplayName', 'Violations Motor 1', 'LineWidth', 2);

    % Improve labels, title, and legend
    xlabel('Time (s)', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Torque (Nm)', 'FontSize', 20, 'FontWeight', 'bold');
    title('Motor 1 Torque Validation', 'FontSize', 22, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 16);
    grid on;
    set(gca, 'FontSize', 24, 'GridColor', tum_grey);

    % Plot results for Motor 2 if applicable
    if config_motor.torque2 > 0
        figure('Position', [100, 100, 1400, 700]);
        plot(time, required_torque_motor2, '-', 'Color', tum_dark_blue, 'LineWidth', 2.5, 'DisplayName', 'Required Torque Motor 2');
        hold on;
        plot(time, available_torque_motor2, '--', 'Color', tum_light_blue, 'LineWidth', 2.5, 'DisplayName', 'Available Torque Motor 2');
        scatter(time(violations_motor2), required_torque_motor2(violations_motor2), 80, 'x', 'MarkerEdgeColor', tum_light_blue, 'DisplayName', 'Violations Motor 2', 'LineWidth', 2);

        % Improve labels, title, and legend
        xlabel('Time (s)', 'FontSize', 20, 'FontWeight', 'bold');
        ylabel('Torque (Nm)', 'FontSize', 20, 'FontWeight', 'bold');
        title('Motor 2 Torque Validation', 'FontSize', 22, 'FontWeight', 'bold');
        legend('Location', 'northeast', 'FontSize', 16);
        grid on;
        set(gca, 'FontSize', 24, 'GridColor', tum_grey);
    end
end