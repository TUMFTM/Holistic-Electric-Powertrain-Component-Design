function [GB, GB_Stufen, config_vehicle, GB_output, config_gearbox] = updateGearboxInformation(config_gearbox, GB_Topology, config_vehicle)
    % UPDATE_GEARBOX_INFORMATION - Updates gearbox properties and integrates them into the vehicle configuration.
    %
    % Inputs:
    %   - config_gearbox: Struct containing gearbox configuration properties.
    %   - GB_Topology: String specifying the gearbox topology ('achsparallel', 'koaxial', or 'einstufig').
    %   - config_vehicle: Struct containing vehicle configuration properties.
    %
    % Outputs:
    %   - GB: Struct containing calculated gearbox parameters.
    %   - GB_Stufen: Number of gearbox stages.
    %   - config_vehicle: Updated vehicle configuration.
    %   - GB_output: Struct containing calculated gearbox losses (gear, seal, bearing, total).
    
    % Initialize gearbox material properties and safety limits
    config_gearbox.power_oilpump = 500; % [W] Maximum power of the oil pump
    config_gearbox.sigma_H_lim = 1500; % [N/mm^2] Fatigue/endurance limit (Naunheimer, 2019, Tab. 7.1)
    config_gearbox.sigma_b_zul = 520 * 0.85 * 0.8 / 1.5; % [N/mm^2] Permissible bending stress (Naunheimer, 2019, p. 478)
    config_gearbox.alpha = 20; % [deg] Normal engagement angle (standard value)
    config_gearbox.beta = [23.4, 18.4]; % [deg] Helix angles (Gao reference)
    config_gearbox.rho_g = 7850; % [kg/m^3] Density of gears (steel)

    % Calculate gearbox design based on topology
    [GB, GB_Stufen] = calculateGearboxDesign(GB_Topology, config_gearbox);

    % Calculate gearbox weight and volume
    % Source: Application of the TOPSIS Method for Multi-Objective Optimization of a Two-Stage Helical Gearbox
    [m_gearbox, V_gearbox] = calcgearboxweight(GB_Topology, GB, config_gearbox);
    GB.m_gearbox = m_gearbox;
    GB.V_gearbox = V_gearbox;

    % Update vehicle configuration with gearbox weight
    config_vehicle = updateGearboxWeight(config_vehicle, GB);

    % Operational parameters for loss calculation
    Mrad = 150; % [Nm] Input torque
    nrad = 3000; % [rpm] Input speed
    Temp = 80; % [°C] Gearbox temperature

    % Calculate gearbox losses based on topology
    switch GB_Topology
        case 'achsparallel'
            [P_V, P_VD, P_VZ, P_VL] = GetrieberechnungFW(GB, Mrad, nrad, Temp);
        case 'einstufig'
            [P_V, P_VD, P_VZ, P_VL] = Getrieberechnung_EinstufigFW(GB, Mrad, nrad, Temp);
        case 'koaxial'
            [P_V, P_VD, P_VZ, P_VL] = Getrieberechnung_koaxialFW(GB, Mrad, nrad, Temp);
        otherwise
            error('Invalid topology. Cannot calculate losses.');
    end

    % Store calculated losses in GB_output
    GB_output.gear_loss = P_VZ;    % Gear loss [W]
    GB_output.seal_loss = P_VD;   % Seal loss [W]
    GB_output.bearing_loss = P_VL; % Bearing loss [W]
    GB_output.total_loss = P_V;   % Total loss [W]
end

function [GB, GB_Stufen] = calculateGearboxDesign(GB_Topology, config_gearbox)
    % CALCULATEGEARBOXDESIGN - Determines the gearbox design based on the topology.
    %
    % Inputs:
    %   - GB_Topology: String specifying the gearbox topology.
    %   - config_gearbox: Struct containing gearbox configuration properties.
    %
    % Outputs:
    %   - GB: Struct containing calculated gearbox parameters.
    %   - GB_Stufen: Number of gearbox stages.

    switch GB_Topology
        case 'achsparallel'
            GB = Getriebedesign(config_gearbox.M_max, config_gearbox.iges, 2); % Parallel-axis
            GB_Stufen = 2;
        case 'koaxial'
            GB = Getriebedesign_koaxial(config_gearbox.M_max, config_gearbox.iges); % Coaxial
            GB_Stufen = 2;
        case 'einstufig'
            GB = Getriebedesign(config_gearbox.M_max, config_gearbox.iges, 1); % Single-stage
            GB_Stufen = 1;
        otherwise
            error('Invalid GB_Topology specified. Choose "achsparallel", "koaxial", or "einstufig".');
    end
end