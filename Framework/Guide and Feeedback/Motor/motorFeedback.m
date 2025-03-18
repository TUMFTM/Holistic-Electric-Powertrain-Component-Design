function motorFeedback(config_motor, avgMotor, config_vehicle, tolerancePercentage, dataTable, secondMotor)
    % MOTOR_USER_FEEDBACK - Provides feedback on motor configuration relative to segment averages and predicted values.
    %
    % Inputs:
    %   config_motor - Structure containing motor configuration details.
    %   avgMotor - Structure containing average motor statistics (e.g., torque, power).
    %   config_vehicle - Structure containing vehicle configuration details.
    %   tolerancePercentage - Tolerance range for comparisons (%).
    %   dataTable - Data table containing vehicle weight and motor statistics for interpolation.
    %
    % Description:
    %   This function compares the user's motor configuration against the average motor 
    %   specifications for a specific segment and predicted values based on vehicle weight.
    %   It highlights any deviations and provides feedback.

    % Calculate the configured motor torque and power
    if secondMotor
        configuredTorque = config_motor.torque1 * config_motor.motorsCount1 + config_motor.torque2 * config_motor.motorsCount2;
        configuredPower = config_motor.power1 * config_motor.motorsCount1 + config_motor.power2 * config_motor.motorsCount2;
    else
        configuredTorque = config_motor.torque1 * config_motor.motorsCount1 ;
        configuredPower = config_motor.power1 * config_motor.motorsCount1 ;
    end    

    % Access average values from avgMotor
    avgValue_MaxTorque = avgMotor.MaxTorque;
    avgValue_MaxPower = avgMotor.MaxPower;

    % Calculate the tolerance range for torque
    lowerLimitTorque = avgValue_MaxTorque * (1 - tolerancePercentage / 100);
    upperLimitTorque = avgValue_MaxTorque * (1 + tolerancePercentage / 100);

    % Display feedback based on the comparison of torque
    if configuredTorque >= lowerLimitTorque && configuredTorque <= upperLimitTorque
        disp(['The configured motor torque of ', num2str(configuredTorque), ...
            ' Nm is close to the average of ', num2str(avgValue_MaxTorque), ...
            ' Nm for the selected segment.']);
    elseif configuredTorque < lowerLimitTorque
        disp(['The configured motor torque of ', num2str(configuredTorque), ...
            ' Nm is significantly lower than the average of ', num2str(avgValue_MaxTorque), ...
            ' Nm for the selected segment. Consider increasing it to match the segment average.']);
    else
        disp(['The configured motor torque of ', num2str(configuredTorque), ...
            ' Nm is significantly higher than the average of ', num2str(avgValue_MaxTorque), ...
            ' Nm for the selected segment. Consider reducing it to avoid excessive performance relative to the segment average.']);
    end

    % Calculate the tolerance range for peak power
    lowerLimitPower = avgValue_MaxPower * (1 - tolerancePercentage / 100);
    upperLimitPower = avgValue_MaxPower * (1 + tolerancePercentage / 100);

    % Display feedback based on the comparison of peak power
    if configuredPower >= lowerLimitPower && configuredPower <= upperLimitPower
        disp(['The configured motor peak power of ', num2str(configuredPower), ...
            ' kW is close to the average of ', num2str(avgValue_MaxPower), ...
            ' kW for the selected segment.']);
    elseif configuredPower < lowerLimitPower
        disp(['The configured motor peak power of ', num2str(configuredPower), ...
            ' kW is significantly lower than the average of ', num2str(avgValue_MaxPower), ...
            ' kW for the selected segment. Consider increasing it to match the segment average.']);
    else
        disp(['The configured motor peak power of ', num2str(configuredPower), ...
            ' kW is significantly higher than the average of ', num2str(avgValue_MaxPower), ...
            ' kW for the selected segment. Consider reducing it to avoid excessive performance relative to the segment average.']);
    end

    % Predict power and torque based on vehicle weight
    vehicleWeight = config_vehicle.mass_total; 
    predictedPower = generalizedInterpolation(dataTable, 'WeightBeforeTeardown_kg_', 'MaxPower_kW_', vehicleWeight);
    predictedTorque = generalizedInterpolation(dataTable, 'WeightBeforeTeardown_kg_', 'MaxTorque_Nm_', vehicleWeight);

    % Define tolerance range for predicted values
    lowerLimitPredictedTorque = predictedTorque * (1 - tolerancePercentage / 100);
    upperLimitPredictedTorque = predictedTorque * (1 + tolerancePercentage / 100);
    lowerLimitPredictedPower = predictedPower * (1 - tolerancePercentage / 100);
    upperLimitPredictedPower = predictedPower * (1 + tolerancePercentage / 100);

    % Display feedback based on the comparison of predicted torque
    if configuredTorque >= lowerLimitPredictedTorque && configuredTorque <= upperLimitPredictedTorque
        disp(['The configured motor torque of ', num2str(configuredTorque), ...
            ' Nm is within the expected range of ', num2str(predictedTorque), ...
            ' Nm based on vehicle weight predictions.']);
    elseif configuredTorque < lowerLimitPredictedTorque
        disp(['The configured motor torque of ', num2str(configuredTorque), ...
            ' Nm is lower than the predicted ', num2str(predictedTorque), ...
            ' Nm for the given vehicle weight. Consider increasing it for expected performance.']);
    else
        disp(['The configured motor torque of ', num2str(configuredTorque), ...
            ' Nm is higher than the predicted ', num2str(predictedTorque), ...
            ' Nm for the given vehicle weight. Ensure this aligns with design requirements.']);
    end

    % Display feedback based on the comparison of predicted power
    if configuredPower >= lowerLimitPredictedPower && configuredPower <= upperLimitPredictedPower
        disp(['The configured motor power of ', num2str(configuredPower), ...
            ' kW is within the expected range of ', num2str(predictedPower), ...
            ' kW based on vehicle weight predictions.']);
    elseif configuredPower < lowerLimitPredictedPower
        disp(['The configured motor power of ', num2str(configuredPower), ...
            ' kW is lower than the predicted ', num2str(predictedPower), ...
            ' kW for the given vehicle weight. Consider increasing it for expected performance.']);
    else
        disp(['The configured motor power of ', num2str(configuredPower), ...
            ' kW is higher than the predicted ', num2str(predictedPower), ...
            ' kW for the given vehicle weight. Ensure this aligns with design requirements.']);
    end
end