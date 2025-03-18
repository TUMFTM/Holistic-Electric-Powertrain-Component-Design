function showPackEquations()
    % Create a figure
    f = figure('Name', 'Battery Pack Equations', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);

    % Add a tab group to organize sections
    tgroup = uitabgroup(f);
    tab1 = uitab(tgroup, 'Title', 'Automatic Design');
    tab2 = uitab(tgroup, 'Title', 'Pack Properties');

    % Default View
    usePlainText = true;

    % Create toggle button
    toggleButton = uicontrol('Style', 'pushbutton', 'String', 'Switch to Parameter View', ...
                             'Position', [650, 550, 120, 30], ...
                             'Callback', @(src, event) toggleView());

    % Axes to add content (text)
    axTab1 = axes('Parent', tab1, 'Position', [0, 0, 1, 1], 'Visible', 'off');
    axTab2 = axes('Parent', tab2, 'Position', [0, 0, 1, 1], 'Visible', 'off');

    % Render initial content
    renderContent(axTab1, axTab2, usePlainText);

    % Function to toggle between views
    function toggleView()
        usePlainText = ~usePlainText;
        if usePlainText
            toggleButton.String = 'Switch to Parameter View';
        else
            toggleButton.String = 'Switch to Text View';
        end
        renderContent(axTab1, axTab2, usePlainText);
    end

    % Function to render content based on view
    function renderContent(axTab1, axTab2, plainText)
        % Clear existing content
        cla(axTab1);
        cla(axTab2);

        if plainText
            % Tab 1: Automatic Design Equations (Plain Text)
            text(axTab1, 0.05, 0.9, "Automatic Design Equations:", 'FontSize', 12);
            text(axTab1, 0.05, 0.8, "Serial Connections:", 'FontSize', 12);
            text(axTab1, 0.1, 0.7, "Number of serial connections = Nominal system voltage ÷ Nominal module voltage", 'FontSize', 10);
            text(axTab1, 0.05, 0.6, "Parallel Connections:", 'FontSize', 12);
            text(axTab1, 0.1, 0.5, "Raw parallel connections = Max(system max current ÷ module max current,", 'FontSize', 10);
            text(axTab1, 0.1, 0.4, "                          system min capacity ÷ module capacity)", 'FontSize', 10);

            % Tab 2: Dependent Pack Properties (Plain Text)
            text(axTab2, 0.05, 0.9, "Dependent Pack Properties:", 'FontSize', 12);
            text(axTab2, 0.05, 0.8, "Electrical Properties:", 'FontSize', 12);
            text(axTab2, 0.1, 0.7, "Pack capacity = Module capacity × Parallel connections", 'FontSize', 10);
            text(axTab2, 0.1, 0.6, "Nominal pack voltage = Nominal module voltage × Serial connections", 'FontSize', 10);
            text(axTab2, 0.1, 0.5, "Pack energy = Pack capacity × Nominal pack voltage ÷ 1,000 (kWh)", 'FontSize', 10);
            text(axTab2, 0.1, 0.4, "Maximum pack voltage = Maximum module voltage × Serial connections", 'FontSize', 10);
            text(axTab2, 0.1, 0.3, "Minimum pack voltage = Minimum module voltage × Serial connections", 'FontSize', 10);
            text(axTab2, 0.05, 0.2, "Dimensional Properties:", 'FontSize', 12);
            text(axTab2, 0.1, 0.15, "Pack mass = Module mass × Total modules × Safety factor for mass", 'FontSize', 10);
            text(axTab2, 0.1, 0.1, "Pack dimensions (x/y/z) = Module dimensions (x/y/z) × Safety factor for dimensions", 'FontSize', 10);
        else
            % Tab 1: Automatic Design Equations (Parameter View)
            text(axTab1, 0.05, 0.9, "Automatic Design Equations:", 'FontSize', 12);
            text(axTab1, 0.05, 0.8, "Serial Connections:", 'FontSize', 12);
            text(axTab1, 0.1, 0.7, "s_sys = ceil(U_sys_nom / U_nom_mod)", 'FontSize', 10);
            text(axTab1, 0.05, 0.6, "Parallel Connections:", 'FontSize', 12);
            text(axTab1, 0.1, 0.5, "p_sys_raw = max(ceil(I_sys_max / I_max_mod), ceil(C_sys_min / C_mod))", 'FontSize', 10);

            % Tab 2: Dependent Pack Properties (Parameter View)
            text(axTab2, 0.05, 0.9, "Dependent Pack Properties:", 'FontSize', 12);
            text(axTab2, 0.05, 0.8, "Electrical Properties:", 'FontSize', 12);
            text(axTab2, 0.1, 0.7, "C_sys = C_mod × p_sys", 'FontSize', 10);
            text(axTab2, 0.1, 0.6, "U_nom_sys = U_nom_mod × s_sys", 'FontSize', 10);
            text(axTab2, 0.1, 0.5, "E_sys = C_sys × U_nom_sys × 10^-3", 'FontSize', 10);
            text(axTab2, 0.1, 0.4, "U_max_sys = U_max_mod × s_sys", 'FontSize', 10);
            text(axTab2, 0.1, 0.3, "U_min_sys = U_min_mod × s_sys", 'FontSize', 10);
            text(axTab2, 0.05, 0.2, "Dimensional Properties:", 'FontSize', 12);
            text(axTab2, 0.1, 0.15, "mass_sys = mass_mod × num_mods_sys × sf_mass_sys", 'FontSize', 10);
            text(axTab2, 0.1, 0.1, "dim_x_sys = dim_x_mod × s_sys × sf_dim_sys", 'FontSize', 10);
            text(axTab2, 0.1, 0.05, "dim_y_sys = dim_y_mod × pe_mod × sf_dim_sys", 'FontSize', 10);
        end
    end
end