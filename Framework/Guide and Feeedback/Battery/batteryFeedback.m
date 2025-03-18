function batteryFeedback(battery_GD, avgBattery,tolerancePercentage)
    % Convert values to numeric if necessary

    avgValue_BatteryEnergy = avgBattery.batteryEnergy;
    avgValue_BatteryCapacity = avgBattery.batteryCapacity;
    avgValue_BatteryVoltage = avgBattery.batteryVoltage;
    avgValue_BatteryWeight = avgBattery.batteryWeight;
    avgValue_BatteryVolume = avgBattery.batteryVolume;
    

    % Calculate battery volume
    configuredVolume = battery_GD.SysInfo.dim_x_sys * battery_GD.SysInfo.dim_y_sys * battery_GD.SysInfo.dim_z_sys;

    % Define tolerance ranges for energy, capacity, weight, and volume
    lowerLimitEnergy = avgValue_BatteryEnergy * (1 - tolerancePercentage / 100);
    upperLimitEnergy = avgValue_BatteryEnergy * (1 + tolerancePercentage / 100);

    lowerLimitCapacity = avgValue_BatteryCapacity * (1 - tolerancePercentage / 100);
    upperLimitCapacity = avgValue_BatteryCapacity * (1 + tolerancePercentage / 100);

    lowerLimitWeight = avgValue_BatteryWeight * (1 - tolerancePercentage / 100);
    upperLimitWeight = avgValue_BatteryWeight * (1 + tolerancePercentage / 100);

    lowerLimitVolume = avgValue_BatteryVolume * (1 - tolerancePercentage / 100);
    upperLimitVolume = avgValue_BatteryVolume * (1 + tolerancePercentage / 100);

    % Feedback for battery energy
    if battery_GD.SysInfo.E_sys >= lowerLimitEnergy && battery_GD.SysInfo.E_sys <= upperLimitEnergy
        disp(['The configured battery energy of ', num2str(battery_GD.SysInfo.E_sys), ...
            ' kWh is close to the average of ', num2str(avgValue_BatteryEnergy), ...
            ' kWh for the selected segment.']);
    elseif battery_GD.SysInfo.E_sys < lowerLimitEnergy
        disp(['The configured battery energy of ', num2str(battery_GD.SysInfo.E_sys), ...
            ' kWh is significantly lower than the average of ', num2str(avgValue_BatteryEnergy), ...
            ' kWh for the selected segment. Consider increasing it to match the segment average.']);
    else
        disp(['The configured battery energy of ', num2str(battery_GD.SysInfo.E_sys), ...
            ' kWh is significantly higher than the average of ', num2str(avgValue_BatteryEnergy), ...
            ' kWh for the selected segment. Consider reducing it to avoid excessive performance relative to the segment average.']);
    end

    % Feedback for battery capacity
    if battery_GD.SysInfo.C_sys >= lowerLimitCapacity && battery_GD.SysInfo.C_sys <= upperLimitCapacity
        disp(['The configured battery capacity of ', num2str(battery_GD.SysInfo.C_sys), ...
            ' Ah is close to the average of ', num2str(avgValue_BatteryCapacity), ...
            ' Ah for the selected segment.']);
    elseif battery_GD.SysInfo.C_sys < lowerLimitCapacity
        disp(['The configured battery capacity of ', num2str(battery_GD.SysInfo.C_sys), ...
            ' Ah is significantly lower than the average of ', num2str(avgValue_BatteryCapacity), ...
            ' Ah for the selected segment. Consider increasing it to match the segment average.']);
    else
        disp(['The configured battery capacity of ', num2str(battery_GD.SysInfo.C_sys), ...
            ' Ah is significantly higher than the average of ', num2str(avgValue_BatteryCapacity), ...
            ' Ah for the selected segment. Consider reducing it to avoid excessive performance relative to the segment average.']);
    end

    % Feedback for battery weight
    if battery_GD.SysInfo.mass_sys >= lowerLimitWeight && battery_GD.SysInfo.mass_sys <= upperLimitWeight
        disp(['The configured battery weight of ', num2str(battery_GD.SysInfo.mass_sys), ...
            ' kg is close to the average of ', num2str(avgValue_BatteryWeight), ...
            ' kg for the selected segment.']);
    elseif battery_GD.SysInfo.mass_sys < lowerLimitWeight
        disp(['The configured battery weight of ', num2str(battery_GD.SysInfo.mass_sys), ...
            ' kg is significantly lower than the average of ', num2str(avgValue_BatteryWeight), ...
            ' kg for the selected segment. Consider increasing it to match the segment average.']);
    else
        disp(['The configured battery weight of ', num2str(battery_GD.SysInfo.mass_sys), ...
            ' kg is significantly higher than the average of ', num2str(avgValue_BatteryWeight), ...
            ' kg for the selected segment. Consider reducing it to avoid excessive performance relative to the segment average.']);
    end

    % Feedback for battery volume
    if configuredVolume >= lowerLimitVolume && configuredVolume <= upperLimitVolume
        disp(['The configured battery volume of ', num2str(configuredVolume), ...
            ' m³ is close to the average of ', num2str(avgValue_BatteryVolume), ...
            ' m³ for the selected segment.']);
    elseif configuredVolume < lowerLimitVolume
        disp(['The configured battery volume of ', num2str(configuredVolume), ...
            ' m³ is significantly lower than the average of ', num2str(avgValue_BatteryVolume), ...
            ' m³ for the selected segment. Consider increasing it to match the segment average.']);
    else
        disp(['The configured battery volume of ', num2str(configuredVolume), ...
            ' m³ is significantly higher than the average of ', num2str(avgValue_BatteryVolume), ...
            ' m³ for the selected segment. Consider reducing it to avoid excessive performance relative to the segment average.']);
    end
end