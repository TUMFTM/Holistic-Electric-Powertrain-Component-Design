function showGearboxEquations()
    % Create the main figure
    f = figure('Name', 'Gearbox Equations', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

    % Add a tab group for the three topologies
    tgroup = uitabgroup(f);
    singleSpeedTab = uitab(tgroup, 'Title', 'Single-Speed Gearbox');
    parallelTab = uitab(tgroup, 'Title', 'Parallel-Axis Gearbox');
    coaxialTab = uitab(tgroup, 'Title', 'Coaxial Gearbox');
    
    % Default View
    usePlainText = true;

    % Create toggle button
    toggleButton = uicontrol('Style', 'pushbutton', 'String', 'Switch to Parameter View', ...
                             'Position', [750, 550, 120, 30], ...
                             'Callback', @(src, event) toggleView());

    % Axes to add content
    axSingleSpeed = axes('Parent', singleSpeedTab, 'Position', [0, 0, 1, 1], 'Visible', 'off');
    axParallel = axes('Parent', parallelTab, 'Position', [0, 0, 1, 1], 'Visible', 'off');
    axCoaxial = axes('Parent', coaxialTab, 'Position', [0, 0, 1, 1], 'Visible', 'off');
    
    % Render initial content
    renderContent(axSingleSpeed, axParallel, axCoaxial, usePlainText);

    % Toggle between plain text and parameterized views
    function toggleView()
        usePlainText = ~usePlainText;
        if usePlainText
            toggleButton.String = 'Switch to Parameter View';
        else
            toggleButton.String = 'Switch to Text View';
        end
        renderContent(axSingleSpeed, axParallel, axCoaxial, usePlainText);
    end

    % Function to render content
    function renderContent(axSingleSpeed, axParallel, axCoaxial, plainText)
        % Clear previous content
        cla(axSingleSpeed);
        cla(axParallel);
        cla(axCoaxial);

        if plainText
            % Single-Speed Gearbox
            text(axSingleSpeed, 0.05, 0.8, "Number of teeth on the wheel = Number of teeth on the pinion × Total ratio", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.75, "Center distance = 0.255 × nthroot(Torque × (Total ratio + 1)^4 / Total ratio, 3)", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.7, "Pitch circle diameter of the pinion = 2 × Center distance / (1 + Number of teeth on the wheel / Number of teeth on the pinion)", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.65, "Pitch circle diameter of the wheel = 2 × Center distance / (1 + Number of teeth on the pinion / Number of teeth on the wheel)", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.4, "Total power loss = Sum of gear losses, bearing losses, and seal losses", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.35, "Gear losses = Load-dependent gear losses + Load-independent gear losses", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.3, "Load-dependent gear losses = Power transmitted × Friction coefficient × Loss factor", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.25, "Load-independent gear losses = Impulse losses + Squeeze losses (if applicable)", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.2, "Bearing losses = Rolling friction moment + Sliding friction moment", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.15, "Seal losses = Calculated using ISO 14179-2", 'FontSize', 10);


            % Parallel-Axis Gearbox
            text(axParallel, 0.05, 0.95, "Stage 1 ratio = 0.7332 × Total ratio^0.6438", 'FontSize', 9);
            text(axParallel, 0.05, 0.9, "Stage 2 ratio = Total ratio / Stage 1 ratio", 'FontSize', 9);
            text(axParallel, 0.05, 0.85, "Teeth on pinion (Stage 1) = Interpolated based on Stage 1 ratio", 'FontSize', 9);
            text(axParallel, 0.05, 0.8, "Teeth on wheel (Stage 1) = Teeth on pinion × Stage 1 ratio", 'FontSize', 9);
            text(axParallel, 0.05, 0.75, "Center distance (Stage 1) = 0.255 × nthroot(Torque × (Stage 1 ratio + 1)^4 / Stage 1 ratio, 3)", 'FontSize', 9);
            text(axParallel, 0.05, 0.7, "Pinion diameter (Stage 1) = 2 × Center distance / (1 + Teeth on wheel / Teeth on pinion)", 'FontSize', 9);
            text(axParallel, 0.05, 0.65, "Wheel diameter (Stage 1) = 2 × Center distance / (1 + Teeth on pinion / Teeth on wheel)", 'FontSize', 9);
            text(axParallel, 0.05, 0.6, "Total power loss = Gear losses + Bearing losses + Seal losses", 'FontSize', 9);
            text(axParallel, 0.05, 0.55, "Gear losses = Load-dependent losses + Load-independent losses", 'FontSize', 9);
            text(axParallel, 0.05, 0.5, "Load-dependent losses = Power × Friction coefficient × Loss factor", 'FontSize', 9);
            text(axParallel, 0.05, 0.45, "Load-independent losses = Impulse losses + Squeeze losses", 'FontSize', 9);
            text(axParallel, 0.05, 0.4, "Bearing losses = Rolling friction moment + Sliding friction moment", 'FontSize', 9);
            text(axParallel, 0.05, 0.35, "Seal losses = Calculated using ISO 14179-2", 'FontSize', 9);
            text(axParallel, 0.05, 0.3, "Gearbox Length = (Pinion Dia. (Stage 1) + Wheel Dia. (Stage 1) / 2 + Pinion Dia. (Stage 2) / 2 + Wheel Dia. (Stage 2) / 2 + 22.5) / 0.975", 'FontSize', 9);
            text(axParallel, 0.05, 0.25, "Gearbox Height = max(Wheel Dia. (Stage 1), Wheel Dia. (Stage 2)) + 8.5 × Support Gap", 'FontSize', 9);
            text(axParallel, 0.05, 0.2, "Gearbox Breadth = Tooth Width (Stage 1) + Tooth Width (Stage 2) + 6 × Support Gap", 'FontSize', 9);
            text(axParallel, 0.05, 0.15, "Gearbox Volume = Length × Height × Breadth", 'FontSize', 9);
            text(axParallel, 0.05, 0.1, "Shafts Volume = Sum of Shaft Volumes", 'FontSize', 9);
            text(axParallel, 0.05, 0.05, "Gearbox Mass = Shafts Mass + Housing Mass + Gears Mass", 'FontSize', 9);

            % Coaxial Gearbox
            text(axCoaxial, 0.05, 0.95, "Stage 1 ratio = 0.299^3 × Torque × (Stage 1 ratio + 1)^4 / Stage 1 ratio = 0.255^3 × Torque × Stage 1 ratio × (Total ratio / Stage 1 ratio + 1)^4 / (Total ratio / Stage 1 ratio)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.9, "Stage 2 ratio = Total ratio / Stage 1 ratio", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.85, "Number of teeth on the pinion for stage 1 = Interpolated based on Stage 1 ratio", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.8, "Number of teeth on the wheel for stage 1 = Number of teeth on the pinion for stage 1 × Stage 1 ratio", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.75, "Hollow shaft outer diameter = Solved symbolically to satisfy bending stress limits", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.7, "Intermediate shaft diameter = 1.72 × nthroot(0.5 × Torque × Total ratio / Allowable shear stress, 3)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.65, "Center distance for stage 1 = 0.299 × nthroot(Torque × (Stage 1 ratio + 1)^4 / Stage 1 ratio, 3)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.6, "Center distance for stage 2 = 0.255 × nthroot(Torque × (Stage 2 ratio + 1)^4 / Stage 2 ratio, 3)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.55, "Tooth width = 4278 × 0.65 × Torque × (Stage ratio + 1) / (Pitch circle diameter^2 × Stage ratio × Contact stress limit^2)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.5, "Normal module = (2 × Center distance × cos(Helix angle)) / (Number of teeth on the pinion + Number of teeth on the wheel)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.45, "Total power loss = Sum of gear losses, bearing losses, and seal losses", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.4, "Gear losses = Load-dependent gear losses + Load-independent gear losses", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.35, "Load-dependent gear losses = Power transmitted × Friction coefficient × Loss factor", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.3, "Load-independent gear losses = Impulse losses + Squeeze losses (if applicable)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.25, "Bearing losses = Rolling friction moment + Sliding friction moment", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.2, "Seal losses = Calculated using ISO 14179-2", 'FontSize', 10);
        else
            % Parameterized equations for Single-Speed
            text(axSingleSpeed, 0.05, 0.8, "z2 = z × i", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.75, "a = 0.255 × nthroot(M × (i + 1)^4 / i, 3)", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.7, "d1 = 2 × a / (1 + z2 / z)", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.65, "d2 = 2 × a / (1 + z / z2)", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.4, "P_loss_total = Σ(P_VZ) + Σ(P_VL) + Σ(P_VD)", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.35, "P_VZ = P_VZP + P_VZ0", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.3, "P_VZP = |P_transmitted| × my × Hv", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.25, "P_VZ0 = P_VQ + P_VI", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.2, "P_VL = (M_rr + M_sl) × ω", 'FontSize', 10);
            text(axSingleSpeed, 0.05, 0.15, "P_VD = 7.69e-6 × (d^2 × n)", 'FontSize', 10);

            % Parameterized equations for Parallel-Axis
            text(axParallel, 0.05, 0.95, "i1 = 0.7332 × iG^0.6438", 'FontSize', 10);
            text(axParallel, 0.05, 0.9, "i2 = iG / i1", 'FontSize', 10);
            text(axParallel, 0.05, 0.85, "z1 = Interpolated value based on i1", 'FontSize', 10);
            text(axParallel, 0.05, 0.8, "z2 = z1 × i1", 'FontSize', 10);
            text(axParallel, 0.05, 0.75, "a = 0.255 × nthroot(M × (i1 + 1)^4 / i1, 3)", 'FontSize', 10);
            text(axParallel, 0.05, 0.7, "d1 = 2 × a / (1 + z2 / z1)", 'FontSize', 10);
            text(axParallel, 0.05, 0.65, "d2 = 2 × a / (1 + z1 / z2)", 'FontSize', 10);
            text(axParallel, 0.05, 0.6, "P_loss_total = Σ(P_VZ) + Σ(P_VL) + Σ(P_VD)", 'FontSize', 10);
            text(axParallel, 0.05, 0.55, "P_VZ = P_VZP + P_VZ0", 'FontSize', 10);
            text(axParallel, 0.05, 0.5, "P_VZP = |P_transmitted| × my × Hv", 'FontSize', 10);
            text(axParallel, 0.05, 0.45, "P_VZ0 = P_VQ + P_VI", 'FontSize', 10);
            text(axParallel, 0.05, 0.4, "P_VL = (M_rr + M_sl) × ω", 'FontSize', 10);
            text(axParallel, 0.05, 0.35, "P_VD = 7.69e-6 × (d^2 × n)", 'FontSize', 10);
            text(axParallel, 0.05, 0.3, "L = (d1_stage1 + d2_stage1 / 2 + d1_stage2 / 2 + d2_stage2 / 2 + 22.5) / 0.975", 'FontSize', 10);
            text(axParallel, 0.05, 0.25, "H = max(d2_stage1, d2_stage2) + 8.5 × S_G", 'FontSize', 10);
            text(axParallel, 0.05, 0.2, "B = b1_stage1 + b1_stage2 + 6 × S_G", 'FontSize', 10);
            text(axParallel, 0.05, 0.15, "V_gearbox = L × H × B", 'FontSize', 10);
            text(axParallel, 0.05, 0.1, "V_shafts = Σ((π / 4) × d_shaft^2 × L_shaft)", 'FontSize', 10);
            text(axParallel, 0.05, 0.05, "m_gearbox = m_shafts + m_housing + m_gears", 'FontSize', 10);

            % Parameterized equations for Coaxial
            % Coaxial Gearbox Equations in Parameterized View
            text(axCoaxial, 0.05, 0.95, "i1 = Solve(0.299^3 × M × (i1 + 1)^4 / i1 = 0.255^3 × M × i1 × (iG / i1 + 1)^4 / (iG / i1))", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.9, "i2 = iG / i1", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.85, "z1 = Interpolated(i1)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.8, "z2 = z1 × i1", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.75, "d1a = Solve(32 × 1.5 × M × 10^3 / (π × σ_b_zul) = (d1a^4 - (dZw + 1)^4) / d1a)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.7, "dZw = 1.72 × nthroot(0.5 × M × 10^3 × iG / 190, 3)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.65, "a1 = 0.299 × nthroot(M × 10^3 × (i1 + 1)^4 / i1, 3)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.6, "a2 = 0.255 × nthroot(M × 10^3 × (i2 + 1)^4 / i2, 3)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.55, "b = 4278 × 0.65 × M × 10^3 × (i + 1) / (d1^2 × i × σ_H_lim^2)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.5, "mn = (2 × a × cos(β)) / (z + z2)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.45, "P_loss_total = Σ(P_VZ) + Σ(P_VL) + Σ(P_VD)", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.4, "P_VZ = P_VZP + P_VZ0", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.35, "P_VZP = |P_transmitted| × my × Hv", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.3, "P_VZ0 = P_VQ + P_VI", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.25, "P_VL = (M_rr + M_sl) × ω", 'FontSize', 10);
            text(axCoaxial, 0.05, 0.2, "P_VD = 7.69e-6 × (d^2 × n)", 'FontSize', 10);
        end
    end
end