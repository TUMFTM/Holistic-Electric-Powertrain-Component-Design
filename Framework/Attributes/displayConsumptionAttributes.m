function displayConsumptionAttributes(simOut, dc, config_vehicle, battery)
    % DISPLAY_CONSUMPTION_ATTRIBUTES - Displays key energy consumption metrics.
    %
    % Outputs:
    %   - Total energy used (kWh)
    %   - Energy consumption per 100 km (kWh/100km)
    %   - Estimated range (km)

    % Calculate driven distance during cycle
    deltaT = 1;
    distance = cumtrapz(dc.speed) * deltaT / 1e3;
    distance = distance(end);

    % Calculate total energy consumption of system (battery discharge + battery loss)
    total_energy_consumption = -sum(simOut.sim_out.battery_power_in.Data, 'omitnan') * deltaT / 1e3 / 3600 ...
        + sum(simOut.sim_out.battery_powerloss.Data, 'omitnan') * deltaT / 1e3 / 3600;
    energy = total_energy_consumption / distance * 100; % kWh/100km

    % Calculate estimated vehicle range
    range = (battery.SysInfo.E_sys / energy) * 100;

    % Display key consumption metrics
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    disp('🔋  **Vehicle Energy Consumption Summary**');
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    fprintf('🔹 **Total Energy Used:** %.2f kWh\n', total_energy_consumption);
    fprintf('🔹 **Energy Consumption per 100 km:** %.2f kWh/100km\n', energy);
    fprintf('🔹 **Estimated Vehicle Range:** %.2f km\n', range);
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
end