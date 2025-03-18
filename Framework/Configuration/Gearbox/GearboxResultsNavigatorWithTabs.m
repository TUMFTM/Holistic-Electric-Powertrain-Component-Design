function GearboxResultsNavigatorWithTabs()
    % Load necessary variables from the base workspace
    if evalin('base', 'exist(''GB'', ''var'') && exist(''config_gearbox'', ''var'') && exist(''GB_Topology'', ''var'')')
        GB = evalin('base', 'GB');
        GB_Topology = evalin('base', 'GB_Topology');
        config_gearbox = evalin('base', 'config_gearbox');
    else
        uialert(uifigure, 'Input variables are not defined. Please run the input function first.', 'Error');
        return;
    end
     
    % Create the main GUI window
    %fig = uifigure('Name', 'Gearbox Results Navigator', 'Position', [500, 100, 800, 600]);
    screenSize = get(0, 'ScreenSize'); % Get the screen size
    figWidth = 1000; % New width
    figHeight = 700; % New height
    figX = (screenSize(3) - figWidth) / 2; % Center X
    figY = (screenSize(4) - figHeight) / 2; % Center Y
    fig = uifigure('Name', 'Gearbox Results Navigator', ...
                   'Position', [figX, figY, figWidth, figHeight]);
    % Add Tab Group
    tabGroup = uitabgroup(fig, 'Position', [20, 20, 760, 560]);

    % Gearbox Properties Tab
    propertiesTab = uitab(tabGroup, 'Title', 'Gearbox Properties');
    displayGearboxProperties(GB, GB_Topology, propertiesTab);

    % Dimensions & Mass Tab
    dimensionsMassTab = uitab(tabGroup, 'Title', 'Dimensions & Mass');
    displayGearboxDimensionsAndMass(GB, GB_Topology, dimensionsMassTab, config_gearbox);

    % Losses Tab
    lossesTab = uitab(tabGroup, 'Title', 'Estimated Losses');
    displayEstimatedLosses(GB, GB_Topology, lossesTab);

    waitfor(fig);
end

function displayGearboxProperties(GB, GB_Topology, parent)
    % Stage 1 Properties
    if strcmp(GB_Topology, 'achsparallel') || strcmp(GB_Topology, 'koaxial')
        z1_stage1 = GB.z1(1); z2_stage1 = GB.z2(1);
        d1_stage1 = GB.d1(1); d2_stage1 = GB.d2(1);
        d_stage1 = GB.d(1); a_stage1 = GB.a(1);
        b_stage1 = GB.b(1); mn_stage1 = GB.mn(1);
    else
        z1_stage1 = GB.z1; z2_stage1 = GB.z2;
        d1_stage1 = GB.d1; d2_stage1 = GB.d2;
        d_stage1 = GB.d; a_stage1 = GB.a;
        b_stage1 = GB.b; mn_stage1 = GB.mn;
    end

    % Stage 2 Properties
    if strcmp(GB_Topology, 'achsparallel') || strcmp(GB_Topology, 'koaxial')
        z1_stage2 = GB.z1(2); z2_stage2 = GB.z2(2);
        d1_stage2 = GB.d1(2); d2_stage2 = GB.d2(2);
        d_stage2 = GB.d(2); a_stage2 = GB.a(2);
        b_stage2 = GB.b(2); mn_stage2 = GB.mn(2);
    else
        z1_stage2 = []; z2_stage2 = [];
        d1_stage2 = []; d2_stage2 = [];
        d_stage2 = []; a_stage2 = [];
        b_stage2 = []; mn_stage2 = [];
    end

    % Stage 1 Display
    uilabel(parent, 'Text', 'Stage 1 Properties:', 'Position', [20, 450, 200, 22], 'FontSize', 14);
    stage1Params = {'Teeth on Pinion (z1):', 'Teeth on Wheel (z2):', ...
                    'Pitch Circle Diameter (d1):', 'Pitch Circle Diameter (d2):', ...
                    'Shaft Diameter (d):', 'Achsabstand (a):', ...
                    'Tooth Width (b):', 'Normalmodul (mn):'};
    stage1Values = {z1_stage1, z2_stage1, d1_stage1, d2_stage1, d_stage1, a_stage1, b_stage1, mn_stage1};
    displayParameters(parent, stage1Params, stage1Values, 420);

    % Stage 2 Display
    uilabel(parent, 'Text', 'Stage 2 Properties:', 'Position', [20, 200, 200, 22], 'FontSize', 14);
    if ~isempty(z1_stage2)
        stage2Params = {'Teeth on Pinion (z1):', 'Teeth on Wheel (z2):', ...
                        'Pitch Circle Diameter (d1):', 'Pitch Circle Diameter (d2):', ...
                        'Shaft Diameter (d):', 'Achsabstand (a):', ...
                        'Tooth Width (b):', 'Normalmodul (mn):'};
        stage2Values = {z1_stage2, z2_stage2, d1_stage2, d2_stage2, d_stage2, a_stage2, b_stage2, mn_stage2};
        displayParameters(parent, stage2Params, stage2Values, 170);
    else
        uilabel(parent, 'Text', 'No second stage for this topology.', ...
            'Position', [20, 170, 400, 22], 'FontSize', 12, 'HorizontalAlignment', 'left');
    end

 
end

function displayGearboxDimensionsAndMass(GB, GB_Topology, parent, config_gearbox)
    % Check if the topology is Parallel-axis
    %if ~strcmp(GB_Topology, 'achsparallel')
    %    uilabel(parent, 'Text', 'Dimensions & Mass calculations are only valid for achsparallel topology.', ...
    %        'Position', [20, 450, 400, 22], 'FontSize', 12, 'HorizontalAlignment', 'left');
    %    return;
    %end
    
    if strcmp(GB_Topology, 'achsparallel')
        % from: Application of the TOPSIS Method for MultiObjective Optimization of a Two-Stage Helical Gearbox
        % Example calculations for demonstration:
        % Gearbox housing dimensions
        d1_stage1 = GB.d1(1); % Diameter of pinion (Stage 1)
        d2_stage1 = GB.d2(1); % Diameter of wheel (Stage 1)
        b1_stage1 = GB.b(1); % Tooth width (Stage 1)
        d1_stage2 = GB.d1(2); % Diameter of pinion (Stage 2)
        d2_stage2 = GB.d2(2); % Diameter of wheel (Stage 2)
        b1_stage2 = GB.b(2); % Tooth width (Stage 2)
        d_shaft_1 = GB.d(1); % Shaft diameter (Stage 1)
        d_shaft_2 = GB.d(2); % Shaft diameter (Stage 2)
        d_shaft_3 = GB.d(3); % Shaft diameter (Stage 3)
    
        % Constants
        rho_g = 7850; % Density of gears (kg/m^3)
        rho_steel = 7850; % Density of steel (kg/m^3)
        rho_alualloy = 2700; % Density of aluminum alloy (kg/m^3)
    
        % Gearbox housing dimensions
        L = (d1_stage1 + d2_stage1 / 2 + d1_stage2 / 2 + d2_stage2 / 2 + 22.5) / 0.975;
        S_G = 0.005 * L + 4.5;
        H = max(d2_stage1, d2_stage2) + 8.5 * S_G;
        B = b1_stage1 + b1_stage2 + 6 * S_G;
         
        % Shaft lengths
        L_1 = B + 1.2 * d_shaft_1;
        L_2 = B;
        L_3 = B + 1.2 * d_shaft_3;
    
        % Volumes of shafts
        V_shaft_1 = (pi / 4) * d_shaft_1^2 * L_1;
        V_shaft_2 = (pi / 4) * d_shaft_2^2 * L_2;
        V_shaft_3 = (pi / 4) * d_shaft_3^2 * L_3;
        V_shafts = (V_shaft_1 + V_shaft_2 + V_shaft_3)/1000000000;
    
        % Volumes of gears
        V_gear_1_stage1 = (pi / 4) * (d1_stage1^2 - d_shaft_1^2) * b1_stage1;
        V_gear_2_stage1 = (pi / 4) * (d2_stage1^2 - d_shaft_2^2) * b1_stage1;
        V_gear_1_stage2 = (pi / 4) * (d1_stage2^2 - d_shaft_2^2) * b1_stage2;
        V_gear_2_stage2 = (pi / 4) * (d2_stage2^2 - d_shaft_3^2) * b1_stage2;
        V_gears = (V_gear_1_stage1 + V_gear_2_stage1 + V_gear_1_stage2 + V_gear_2_stage2)/1000000000;
    
        % Housing volumes
        V_gearbox = L * H * B/1000000000;
        V_A = L * H * S_G/1000000000;
        V_B = L * b1_stage1 * 1.5 * S_G/1000000000;
        V_C = (B - 2 * S_G) * H * S_G/1000000000;
    
        % Mass calculations
        e_1 = 1; e_2 = 0.6;
        m_g1 = rho_g * ((pi * e_1 * d1_stage1^2 * b1_stage1) / 4 + (pi * e_2 * d2_stage1^2 * b1_stage1) / 4)/1000000000;
        m_g2 = rho_g * ((pi * e_1 * d1_stage2^2 * b1_stage2) / 4 + (pi * e_2 * d2_stage2^2 * b1_stage2) / 4)/1000000000;
        m_shafts = rho_steel * V_shafts;
        m_housing = rho_alualloy * 2 * (V_A + V_B + V_C);
        m_gears = m_g1 + m_g2;
        m_gearbox = m_shafts + m_housing + m_gears;
        
        % Dimensions & Mass Display
        params = {'Shaft Mass (kg):', 'Housing Mass (kg):', 'Gears Mass (kg):', 'Total Mass (kg):', ...
                  'Shaft Volume (m³):', 'Total Volume (m³):'};
        values = {m_shafts, m_housing, m_gears, m_gearbox, V_shafts, V_gearbox};
    elseif strcmp(GB_Topology, 'koaxial')
        m_gearbox = 0.199*((GB.iges*config_gearbox.M_max)^0.669*GB.Stufen^0.334);
        V_gearbox = 0;
        params = {'Total Mass (kg):'};
        values = {m_gearbox};
    else
        m_gearbox = 0.08*config_gearbox.M_max^0.75*GB.iges^0.2;
        V_gearbox =0;
        params = {'Total Mass (kg):'};
        values = {m_gearbox};
    end
        
        GB.m_gearbox = m_gearbox;
        %GB_output.V_gearbox = V_gearbox;
        assignin('base', 'GB', GB);
        %assignin('base', 'GB_output', GB_output);
        
        displayParameters(parent, params, values, 450);
end

function displayEstimatedLosses(GB, GB_Topology, parent)
     % Fixed input values for loss estimation
    Mrad = 150; % Nm
    nrad = 3000; % rpm
    Temp = 80; % °C

    %GB_output = evlain('base', 'GB_output');
    switch GB_Topology
        case 'achsparallel'
            [P_V, P_VD, P_VZ, P_VL] = GetrieberechnungFW(GB, Mrad, nrad, Temp);
        case 'einstufig'
            [P_V, P_VD, P_VZ, P_VL] = Getrieberechnung_EinstufigFW(GB, Mrad, nrad, Temp);
        case 'koaxial'
            [P_V, P_VD, P_VZ, P_VL] = Getrieberechnung_koaxialFW(GB, Mrad, nrad, Temp);
        otherwise
            uialert([], 'Invalid topology. Cannot calculate losses.', 'Error');
            return;
    end

    % Losses Display
    params = {'Verzahnungsverluste Leistung (P_{VZ}):', ...
              'Dichtungsverluste Leistung (P_{VD}):', ...
              'Lagerverluste Leistung (P_{VL}):', ...
              'Gesamtverluste Leistung (P_{V}):'};
    values = {P_VZ, P_VD, P_VL, P_V};
    
    %GB_output.gear_loss = P_VZ;
    %GB_output.seal_loss = P_VD;
    %GB_output.bearing_loss = P_VL;
    %GB_output.total_loss = P_V;
    %assignin('base', 'GB_output', GB_output);

    displayParameters(parent, params, values, 450);
end

function displayParameters(parent, params, values, yStart)
    col1X = 20; % Left column X
    col2X = 360; % Right column X
    colWidth = 320; % Set a consistent width for each column
    for i = 1:length(params)
        if mod(i, 2) == 1 % Left column
            xPos = col1X;
            yPos = yStart - floor((i - 1) / 2) * 30;
        else % Right column
            xPos = col2X;
            yPos = yStart - floor((i - 1) / 2) * 30;
        end
        uilabel(parent, 'Text', params{i}, 'Position', [xPos, yPos, colWidth - 100, 22], ...
            'FontSize', 12, 'HorizontalAlignment', 'right');
        uilabel(parent, 'Text', num2str(values{i}), 'Position', [xPos + colWidth - 90, yPos, 100, 22], ...
            'FontSize', 12, 'HorizontalAlignment', 'left');
    end
end