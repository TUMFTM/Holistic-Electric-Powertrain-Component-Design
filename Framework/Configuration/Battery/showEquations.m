function showEquations()
    % Create a figure
    f = figure('Name', 'Battery Module Equations', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
    
    % Add a tab group to organize sections
    tgroup = uitabgroup(f);
    tab1 = uitab(tgroup, 'Title', 'Automatic Design');
    tab2 = uitab(tgroup, 'Title', 'Module Properties');
    
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
            text(axTab1, 0.1, 0.7, "Raw serial connections = Min(max serial connections, nominal module voltage ÷ max cell voltage)", 'FontSize', 10);
            text(axTab1, 0.1, 0.6, "Serial connections = Raw serial connections rounded to the nearest even number", 'FontSize', 10);
            text(axTab1, 0.05, 0.5, "Parallel Connections:", 'FontSize', 12);
            text(axTab1, 0.1, 0.4, "Min parallel connections = Max(max system current ÷ max cell current × safety factor,", 'FontSize', 10);
            text(axTab1, 0.1, 0.3, "                          min module capacity ÷ cell capacity)", 'FontSize', 10);
            
            % Tab 2: Dependent Module Properties (Plain Text)
            text(axTab2, 0.05, 0.9, "Dependent Module Properties:", 'FontSize', 12);
            text(axTab2, 0.05, 0.8, "Electrical Properties:", 'FontSize', 12);
            text(axTab2, 0.1, 0.7, "Module capacity = Cell capacity × Parallel connections", 'FontSize', 10);
            text(axTab2, 0.1, 0.6, "Nominal module voltage = Nominal cell voltage × Serial connections", 'FontSize', 10);
            text(axTab2, 0.1, 0.5, "Module energy = Module capacity × Nominal module voltage ÷ 1,000 (kWh)", 'FontSize', 10);
            text(axTab2, 0.1, 0.4, "Maximum module voltage = Maximum cell voltage × Serial connections", 'FontSize', 10);
            text(axTab2, 0.1, 0.3, "Minimum module voltage = Minimum cell voltage × Serial connections", 'FontSize', 10);
            text(axTab2, 0.1, 0.2, "Maximum module current = Maximum cell current × Parallel connections", 'FontSize', 10);
        else
            % Tab 1: Automatic Design Equations (Parameter View)
            text(axTab1, 0.05, 0.9, "Automatic Design Equations:", 'FontSize', 12);
            text(axTab1, 0.05, 0.8, "Serial Connections:", 'FontSize', 12);
            text(axTab1, 0.1, 0.7, "s_mod_raw = min(s_mod_max, floor(U_mod_nom / U_max_cell))", 'FontSize', 10);
            text(axTab1, 0.1, 0.6, "s_mod = 2 × floor(s_mod_raw / 2)", 'FontSize', 10);
            text(axTab1, 0.05, 0.5, "Parallel Connections:", 'FontSize', 12);
            text(axTab1, 0.1, 0.4, "p_min_mod = max(ceil(I_mod_max × (1 + sf_parallel) / I_max_cell), ceil(C_mod_min / C_A_cell))", 'FontSize', 10);
            
            % Tab 2: Dependent Module Properties (Parameter View)
            text(axTab2, 0.05, 0.9, "Dependent Module Properties:", 'FontSize', 12);
            text(axTab2, 0.05, 0.8, "Electrical Properties:", 'FontSize', 12);
            text(axTab2, 0.1, 0.7, "C_mod = C_A_cell × p_mod", 'FontSize', 10);
            text(axTab2, 0.1, 0.6, "U_nom_mod = U_nom_cell × s_mod", 'FontSize', 10);
            text(axTab2, 0.1, 0.5, "E_mod = C_mod × U_nom_mod × 10^-3", 'FontSize', 10);
            text(axTab2, 0.1, 0.4, "U_max_mod = U_max_cell × s_mod", 'FontSize', 10);
            text(axTab2, 0.1, 0.3, "U_min_mod = U_min_cell × s_mod", 'FontSize', 10);
            text(axTab2, 0.1, 0.2, "I_max_mod = I_max_cell × p_mod", 'FontSize', 10);
        end
    end
end