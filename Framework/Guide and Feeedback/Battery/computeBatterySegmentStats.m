function [avgBattery] = computeBatterySegmentStats(segment, dataTable)
    % COMPUTEBATTERYSEGMENTSTATS - Computes average battery statistics for a given vehicle segment.
    %
    % Outputs:
    %   - Console feedback with segment-based battery statistics.

    % Initialize output structure
    avgBattery = struct();

    % Calculate averages for different battery parameters
    [~, avgValue_batteryVoltage, num_batteryVoltage] = calculateSumAndAverage(dataTable, segment, 'BatteryPackNominalVoltage_V__Overview_'); 
    avgBattery.batteryVoltage = avgValue_batteryVoltage;

    [~, avgValue_batteryCapacity, num_batteryCapacity] = calculateSumAndAverage(dataTable, segment, 'BatteryPackCapacity_Ah__Overview_'); 
    avgBattery.batteryCapacity = avgValue_batteryCapacity;

    [~, avgValue_batteryEnergy, num_batteryEnergy] = calculateSumAndAverage(dataTable, segment, 'BatteryPackEnergy_kWh__Overview_'); 
    avgBattery.batteryEnergy = avgValue_batteryEnergy;

    [~, avgValue_batteryWeight, num_batteryWeight] = calculateSumAndAverage(dataTable, segment, 'BatteryPackWeight_Kg___Overview_'); 
    avgBattery.batteryWeight = avgValue_batteryWeight;

    [~, avgValue_batteryVolume, num_batteryVolume] = calculateSumAndAverage(dataTable, segment, 'BatteryPack3DVolume_L__Overview_'); 
    avgValue_batteryVolume = avgValue_batteryVolume / 1000; % Convert from liters to cubic meters
    avgBattery.batteryVolume = avgValue_batteryVolume;

    % Display battery segment statistics
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    fprintf('🔋  **Battery Statistics for Segment %s**\n', segment);
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    fprintf('⚡ **Avg. Nominal Voltage:** %.2f V (from %d vehicles)\n', avgValue_batteryVoltage, num_batteryVoltage);
    fprintf('🔄 **Avg. Capacity:** %.2f Ah (from %d vehicles)\n', avgValue_batteryCapacity, num_batteryCapacity);
    fprintf('⚡ **Avg. Energy Content:** %.2f kWh (from %d vehicles)\n', avgValue_batteryEnergy, num_batteryEnergy);
    fprintf('⚖️ **Avg. Battery Weight:** %.2f kg (from %d vehicles)\n', avgValue_batteryWeight, num_batteryWeight);
    fprintf('📦 **Avg. Battery Volume:** %.4f m³ (from %d vehicles)\n', avgValue_batteryVolume, num_batteryVolume);
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
end