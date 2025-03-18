function validateBattery(simOut, battery_GD)
    % VALIDATE_BATTERY - Validates whether the battery is sufficient for the power and current demands of the motor.
    %
    % Inputs:
    %   - simOut: Simulink simulation output containing motor data.
    %   - battery_GD: Struct containing battery parameters (voltage, current, etc.).
    %
    % Outputs:
    %   - Console feedback with detailed analysis of battery performance.
    %   - Visualization of motor power demand vs. battery capacity.

    % Extract motor power demand from simulation logs
    motor_power = simOut.sim_out.terminal_power.Data;  % Power in Watts
    time = simOut.sim_out.Motor_rpm.Time;  % Time vector for reference

    % Extract battery parameters
    V_battery_nominal = battery_GD.SysInfo.U_nom_sys;  % Nominal voltage (V)
    I_battery_max = battery_GD.SysInfo.I_max_sys;  % Maximum allowable current (A)
    P_battery_max = V_battery_nominal * I_battery_max;  % Maximum power output (W)

    % Compute key statistics for motor power demand
    P_max_required = max(motor_power);  % Max power demand (W)
    P_avg_required = mean(motor_power);  % Average power demand (W)
    P_spare = P_battery_max - P_max_required;  % Spare power capacity (W)

    % Compute the required current at max power
    I_required_max = P_max_required / V_battery_nominal;  % Max required current (A)
    I_spare = I_battery_max - I_required_max;  % Spare current capacity (A)

    % Determine power sufficiency message
    if P_max_required > P_battery_max
        power_message = sprintf('⚠️  Power NOT sufficient!  Required power exceeds battery limit by **%.2f kW**.', ...
                                (P_max_required - P_battery_max) / 1000);
    else
        power_message = sprintf('✅  Power sufficient!  Battery has **%.2f kW** spare capacity.', P_spare / 1000);
    end

    % Determine current sufficiency message
    if I_required_max > I_battery_max
        current_message = sprintf('⚠️  Current NOT sufficient!  Required current exceeds battery limit by **%.2f A**.', ...
                                  I_required_max - I_battery_max);
    else
        current_message = sprintf('✅  Current sufficient!  Battery has **%.2f A** spare capacity.', I_spare);
    end

    % Display battery validation results with cleaner formatting
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    disp('🔋  **Battery Validation Results**');
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    fprintf('🔹 **Max Required Power:**      %.2f kW\n', P_max_required / 1000);
    fprintf('🔹 **Battery Max Power:**       %.2f kW\n', P_battery_max / 1000);
    fprintf('🔹 **Max Required Current:**    %.2f A\n', I_required_max);
    fprintf('🔹 **Battery Max Current:**     %.2f A\n', I_battery_max);
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    disp(power_message);
    disp(current_message);
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');

    % Define TUM color scheme
    tum_blue = [0, 101, 189] / 255;
    tum_dark_blue = [0, 79, 158] / 255;
    tum_light_blue = [100, 160, 200] / 255;
    tum_grey = [162, 162, 162] / 255;

    % Plot motor power demand vs. battery capacity
    figure('Name', 'Battery Operation Validation', 'Position', [100, 100, 1400, 700]);
    plot(time, motor_power / 1000, '-', 'Color', tum_dark_blue, 'LineWidth', 2.5, 'DisplayName', 'Motor Power Demand (kW)');
    hold on;
    yline(P_battery_max / 1000, '--', 'Color', tum_light_blue, 'LineWidth', 2.5, 'DisplayName', 'Battery Maximum Power (kW)');

    % Improve labels, title, and legend
    xlabel('Time (s)', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Power (kW)', 'FontSize', 20, 'FontWeight', 'bold');
    title('Battery Power Demand vs. Capacity', 'FontSize', 22, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 16);
    grid on;

    % Ensure axis font sizes update correctly
    set(gca, 'FontSize', 24, 'FontWeight', 'bold', 'GridColor', tum_grey);
end