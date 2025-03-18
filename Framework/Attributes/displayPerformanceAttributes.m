function max_velocity =displayPerformanceAttributes(simOut, config_gearbox, config_vehicle, motormap, secondMotor)
    % DISPLAY_PERFORMANCE_ATTRIBUTES - Displays acceleration and speed performance results.
    %
    % Outputs:
    %   - Console feedback with acceleration and speed details.
    %   - Graphical plot of velocity over time.

    % Extract time and velocity data from simulation
    time = simOut.sim_out.time.Data;
    velocity = simOut.sim_out.velocity.Data .* 3.6;  % Convert to km/h
    index = find(velocity > 100, 1);
    acceleration_time = time(index);

    % Calculate Torque-Based Peak Speed
    torque_based_peak_speed = velocity(end);

    % Retrieve maximum motor RPM for one or two motors
    max_motor_rpm1 = max(max(motormap.Speed)); % Highest available RPM for motor 1
    if secondMotor
        motormap2 = evalin("base", "motormap2");
        max_motor_rpm2 = max(max(motormap2.Speed)); % Highest available RPM for motor 2
        max_motor_rpm = min(max_motor_rpm1, max_motor_rpm2); % Use the lower RPM limit
    else
        max_motor_rpm = max_motor_rpm1;
    end

    % Retrieve transmission ratio and tire radius
    i_gear = config_gearbox.iges; % Transmission ratio
    r_tire = config_vehicle.tire_radius; % Tire radius (m)

    % Calculate RPM-Based Maximum Velocity
    rpm_based_peak_speed = (max_motor_rpm / i_gear) * (2 * pi / 60) * r_tire * 3.6;
    max_velocity = min(rpm_based_peak_speed, torque_based_peak_speed);
    % Display performance attributes
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    disp('🚀  **Vehicle Performance Summary**');
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    fprintf('🔹 **Acceleration Time (0-100 km/h):**   %.2f seconds\n', acceleration_time);
    fprintf('🔹 **Torque-Based Peak Speed:**          %.2f km/h\n', torque_based_peak_speed);
    fprintf('🔹 **RPM-Based Peak Speed:**             %.2f km/h\n', rpm_based_peak_speed);
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');

    % Define TUM colors
    tum_dark_blue = [0, 79, 158] / 255;     % TUM Dark Blue
    tum_grey = [162, 162, 162] / 255;       % TUM Grey

    % Plot velocity over time
    figure('Name', 'Performance', 'Position', [100, 100, 1400, 700]);
    plot(time, velocity, '-', 'Color', tum_dark_blue, 'LineWidth', 2.5, 'DisplayName', 'Velocity');

    % Improve labels, title, and grid settings
    xlabel('Time (s)', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Velocity (km/h)', 'FontSize', 20, 'FontWeight', 'bold');
    title('Acceleration Cycle', 'FontSize', 22, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 24, 'FontWeight', 'bold', 'GridColor', tum_grey);
end