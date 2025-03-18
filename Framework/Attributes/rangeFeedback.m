function rangeFeedback(dataTable, segment, battery, comparisonFactor, user_range)
    % RANGE_FEEDBACK - Compares vehicle range with segment averages and battery capacity predictions.
    %
    % Outputs:
    %   - Console feedback on range performance.

    % Validate comparisonFactor input
    if comparisonFactor <= 0 || comparisonFactor >= 1
        error('comparisonFactor must be between 0 and 1.');
    end

    upperLimit = 1 + comparisonFactor;
    lowerLimit = 1 - comparisonFactor;

    % Get segment average range from the databank
    [~, avgValue_Range, num_Range] = calculateSumAndAverage(dataTable, segment, 'ElectricDrivingRange_km_OEMData_');

    % Display range comparison summary
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    disp('🔋  **Range Performance Feedback**');
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    fprintf('🔹 **User Vehicle Range:** %.2f km\n', user_range);
    fprintf('🔹 **Segment Average Range:** %.2f km (Based on %d vehicles)\n', avgValue_Range, num_Range);

    if user_range > avgValue_Range * upperLimit
        disp('🔥 **Exceptional Range!** Your vehicle significantly outperforms the segment average.');
    elseif user_range < avgValue_Range * lowerLimit
        disp('⚠️  **Limited Range!** Your vehicle has a shorter range than the segment average. Consider optimizing efficiency or battery capacity.');
    else
        disp('✅ **Balanced Performance!** Your vehicle’s range aligns well with segment expectations.');
    end

    % Cross-check range against battery energy capacity
    calculated_range = generalizedInterpolation(dataTable, 'BatteryPackEnergyOEM_kWh__Overview_', 'ElectricDrivingRange_km_OEMData_', battery.SysInfo.E_sys);

    fprintf('\n🔹 **Predicted Range Based on Battery Capacity (%.2f kWh):** %.2f km\n', battery.SysInfo.E_sys, calculated_range);
    fprintf('🔹 **Simulated User Range:** %.2f km\n', user_range);

    if user_range > calculated_range * upperLimit
        disp('🔥 **Above Expectations!** Your vehicle achieves a longer range than predicted based on battery energy capacity.');
    elseif user_range < calculated_range * lowerLimit
        disp('⚠️  **Underperforming!** Your vehicle’s range is lower than expected. Consider improving efficiency.');
    else
        disp('✅ **Well-Calibrated!** Your vehicle’s range is in line with its battery energy capacity.');
    end

    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
end