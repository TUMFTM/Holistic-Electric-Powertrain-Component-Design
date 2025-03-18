function performanceFeedback(dataTable, segment, simOut, config_motor, comparisonFactor, secondMotor, user_max_velocity)
    % PERFORMANCE_FEEDBACK - Compares vehicle speed and acceleration with segment averages
    %
    % Outputs:
    %   - Console feedback with performance analysis.

    % Validate comparisonFactor input
    if comparisonFactor <= 0 || comparisonFactor >= 1
        error('comparisonFactor must be between 0 and 1.');
    end

    upperLimit = 1 + comparisonFactor;
    lowerLimit = 1 - comparisonFactor;

   

    % Get average top speed from the segment
    [~, avgValue_TopSpeed, num_TopSpeed] = calculateSumAndAverage(dataTable, segment, 'TopSpeed_km_h_');

    % Display performance summary
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    disp('🚀  **Performance Feedback**');
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');

    % Display top speed comparison
    fprintf('🔹 **Maximum Speed:** %.2f km/h\n', user_max_velocity);
    fprintf('🔹 **Segment Average Speed:** %.2f km/h (Based on %d vehicles)\n', avgValue_TopSpeed, num_TopSpeed);

    if user_max_velocity > avgValue_TopSpeed * upperLimit
        disp('🔥 Your vehicle is a **beast!** It exceeds the segment average and leaves competitors behind!');
    elseif user_max_velocity < avgValue_TopSpeed * lowerLimit
        disp('⚠️  Your vehicle lags behind the segment average. It might struggle to keep up.');
    else
        disp('✅ Your vehicle is well-balanced, performing within the expected range for this segment.');
    end

    % Extract acceleration data
    time = simOut.sim_out.time.Data;
    velocity = simOut.sim_out.velocity.Data .* 3.6;
    index = find(velocity > 100, 1);
    user_acceleration_time = time(index);

    % Get average acceleration from the segment
    [~, avgValue_Acceleration, num_Acceleration] = calculateSumAndAverage(dataTable, segment, 'Acceleration0_100Km_h_s_');

    % Display acceleration comparison
    fprintf('\n🔹 **Acceleration (0-100 km/h):** %.2f seconds\n', user_acceleration_time);
    fprintf('🔹 **Segment Average Acceleration:** %.2f seconds (Based on %d vehicles)\n', avgValue_Acceleration, num_Acceleration);

    if user_acceleration_time < avgValue_Acceleration * lowerLimit
        disp('⚡ **Lightning fast!** Your vehicle crushes the segment average!');
    elseif user_acceleration_time > avgValue_Acceleration * upperLimit
        disp('🐢 Your vehicle takes its time. Consider performance tuning.');
    else
        disp('✅ Your vehicle accelerates within the expected range for this segment.');
    end

    % Determine total motor power and torque
    if secondMotor
        motors_total_power = config_motor.power1 * config_motor.motorsCount1 + config_motor.power2 * config_motor.motorsCount2;
        motors_total_torque = config_motor.torque1 * config_motor.motorsCount1 + config_motor.torque2 * config_motor.motorsCount2;
    else
        motors_total_power = config_motor.power1 * config_motor.motorsCount1;
        motors_total_torque = config_motor.torque1 * config_motor.motorsCount1;
    end

    % Predicted vs. actual top speed (based on installed motor power)
    predictedSpeed = generalizedInterpolation(dataTable, 'MaxPower_kW_', 'TopSpeed_km_h_', motors_total_power);
    fprintf('\n🔹 **Predicted Top Speed (Based on Installed Motor Power: %.2f kW):** %.2f km/h\n', motors_total_power, predictedSpeed);
    fprintf('🔹 **Actual Top Speed:** %.2f km/h\n', user_max_velocity);

    if user_max_velocity > predictedSpeed * upperLimit
        disp('🔥 Your vehicle **exceeds expectations!** It outperforms the predicted top speed!');
    elseif user_max_velocity < predictedSpeed * lowerLimit
        disp('⚠️  Your vehicle **falls short** of the predicted top speed. Consider revising your design.');
    else
        disp('✅ Your vehicle performs **as expected**, achieving a speed close to predictions.');
    end

    % Predicted vs. actual acceleration time (based on installed motor torque)
    predicted_acceleration_time = generalizedInterpolation(dataTable, 'MaxTorque_Nm_', 'Acceleration0_100Km_h_s_', motors_total_torque);
    fprintf('\n🔹 **Predicted Acceleration (Based on Installed Motor Torque: %.2f Nm):** %.2f seconds\n', motors_total_torque, predicted_acceleration_time);
    fprintf('🔹 **Actual Acceleration:** %.2f seconds\n', user_acceleration_time);

    if user_acceleration_time < predicted_acceleration_time * lowerLimit
        disp('⚡ Your vehicle **accelerates like a rocket!** Faster than expected!');
    elseif user_acceleration_time > predicted_acceleration_time * upperLimit
        disp('🐢 Your vehicle **struggles** to meet expectations. Optimization may be needed.');
    else
        disp('✅ Your vehicle’s acceleration aligns with predictions. A **finely tuned machine!**');
    end

    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
end