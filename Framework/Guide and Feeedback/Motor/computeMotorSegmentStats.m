function [avgMotor, normPower] = computeMotorSegmentStats(segment, dataTable)
    % COMPUTEMOTORSEGMENTSTATS - Computes motor-related statistics for a given vehicle segment.
    %
    % Outputs:
    %   - Console feedback with motor stats and suggested power normalization.

    % Initialize output structures
    avgMotor = struct();
    normPower = struct();

    %% --- Compute Averages for Motor Parameters ---
    
    % Average front motor power
    [~, avgFrontPower, numFrontPower] = calculateSumAndAverage(dataTable, segment, 'FrontE_Motor1Power_kW_');
    avgMotor.FrontPower = avgFrontPower;

    % Average rear motor power
    [~, avgRearPower, numRearPower] = calculateSumAndAverage(dataTable, segment, 'RearE_Motor1Power_kW_');
    avgMotor.RearPower = avgRearPower;

    % Average maximum power
    [~, avgMaxPower, numMaxPower] = calculateSumAndAverage(dataTable, segment, 'MaxPower_kW_');
    avgMotor.MaxPower = avgMaxPower;

    % Average maximum torque
    [~, avgMaxTorque, numMaxTorque] = calculateSumAndAverage(dataTable, segment, 'MaxTorque_Nm_');
    avgMotor.MaxTorque = avgMaxTorque;

    % Average front motor torque
    [~, avgFrontTorque, numFrontTorque] = calculateSumAndAverage(dataTable, segment, 'FrontE_Motor1Torque_Nm_');
    avgMotor.FrontTorque = avgFrontTorque;

    % Average rear motor torque
    [~, avgRearTorque, numRearTorque] = calculateSumAndAverage(dataTable, segment, 'RearE_Motor1Torque_Nm_');
    avgMotor.RearTorque = avgRearTorque;

    %% --- Count Motor Configurations ---
    motorCountDistribution = countValuesInSegment(dataTable, segment, 'NumberOfE_motor', {1, 2, 3, 4});
    avgMotor.numMotors_1 = motorCountDistribution.x1;
    avgMotor.numMotors_2 = motorCountDistribution.x2;
    avgMotor.numMotors_3 = motorCountDistribution.x3;
    avgMotor.numMotors_4 = motorCountDistribution.x4;

    [~, mostCommonMotorsNum] = max([motorCountDistribution.x1, motorCountDistribution.x2, ...
                                    motorCountDistribution.x3, motorCountDistribution.x4]);
    avgMotor.typicalMotorsNum = mostCommonMotorsNum; 

    % Count front & rear motor types
    motorTypeDistributionFront = countValuesInSegment(dataTable, segment, 'FrontE_Motor1Type', ...
        {'Permanent Magnet / AC Synchronous', 'Induction / AC Asynchronous'});
    motorTypeDistributionRear = countValuesInSegment(dataTable, segment, 'RearE_Motor1Type', ...
        {'Permanent Magnet/AC Synchronous', 'Induction / AC Asynchronous'});

    %% --- Compute Normalized Power Suggestion ---
    if (avgFrontPower + avgRearPower) > 0
        normFrontPower = (avgFrontPower / (avgFrontPower + avgRearPower)) * avgMaxPower;
        normRearPower = (avgRearPower / (avgFrontPower + avgRearPower)) * avgMaxPower;
    else
        normFrontPower = avgFrontPower;
        normRearPower = avgRearPower;
    end

    if (avgFrontTorque + avgRearTorque) > 0
        normFrontTorque = (avgFrontTorque / (avgFrontTorque + avgRearTorque)) * avgMaxTorque;
        normRearTorque = (avgRearTorque / (avgFrontTorque + avgRearTorque)) * avgMaxTorque;
    else
        normFrontTorque = avgFrontTorque;
        normRearTorque = avgRearTorque;
    end

    normPower.FrontPower = normFrontPower;
    normPower.RearPower = normRearPower;
    normPower.TotalPower = avgMaxPower;
    
    normTorque.FrontTorque = normFrontTorque;
    normTorque.RearTorque = normRearTorque;
    normTorque.TotalTorque = avgMaxTorque;

    %% --- Display Results ---
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    fprintf('🚗  **Motor Statistics for Segment %s**\n', segment);
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');

    % Display motor power statistics
    fprintf('🔹 **Avg. Front Motor Power:** %.2f kW (%d vehicles)\n', avgFrontPower, numFrontPower);
    fprintf('🔹 **Avg. Rear Motor Power:** %.2f kW (%d vehicles)\n', avgRearPower, numRearPower);
    fprintf('🔹 **Avg. Max Power:** %.2f kW (%d vehicles)\n', avgMaxPower, numMaxPower);
    fprintf('🔹 **Avg. Max Torque:** %.2f Nm (%d vehicles)\n', avgMaxTorque, numMaxTorque);

    % Display torque statistics
    fprintf('\n🔹 **Avg. Front Motor Torque:** %.2f Nm (%d vehicles)\n', avgFrontTorque, numFrontTorque);
    fprintf('🔹 **Avg. Rear Motor Torque:** %.2f Nm (%d vehicles)\n', avgRearTorque, numRearTorque);

    % Display number of motors distribution
    fprintf('\n🔹 **Motor Configuration Distribution:**\n');
    fprintf('   - %d vehicles have **1 motor**\n', motorCountDistribution.x1);
    fprintf('   - %d vehicles have **2 motors**\n', motorCountDistribution.x2);
    fprintf('   - %d vehicles have **3 motors**\n', motorCountDistribution.x3);
    fprintf('   - %d vehicles have **4 motors**\n', motorCountDistribution.x4);

    % Display front motor type distribution
    disp('🔹 **Front Motor Types:**');
    fieldsFront = fieldnames(motorTypeDistributionFront);
    for i = 1:length(fieldsFront)
        fprintf('   - %d vehicles use "%s"\n', motorTypeDistributionFront.(fieldsFront{i}), fieldsFront{i});
    end

    % Display rear motor type distribution
    disp('🔹 **Rear Motor Types:**');
    fieldsRear = fieldnames(motorTypeDistributionRear);
    for i = 1:length(fieldsRear)
        fprintf('   - %d vehicles use "%s"\n', motorTypeDistributionRear.(fieldsRear{i}), fieldsRear{i});
    end

    % Display normalized power suggestions
    disp('⚡ **Suggested Normalized Power & Torque**');
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    fprintf('🔹 **Front Motor Power:** %.2f kW\n', normFrontPower);
    fprintf('🔹 **Rear Motor Power:** %.2f kW\n', normRearPower);
    fprintf('🔹 **Front Motor Torque:** %.2f Nm\n', normFrontTorque);
    fprintf('🔹 **Rear Motor Torque:** %.2f Nm\n', normRearTorque);
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
end