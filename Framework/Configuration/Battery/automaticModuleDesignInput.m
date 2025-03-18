function [battery_GD, SysSpec] = automaticModuleDesignInput(battery_GD, SysSpec)
    % automaticDesignCallback
    % Creates a GUI to accept user inputs for automatic design parameters
    % and processes battery module configurations based on the provided inputs.
    %
    % Inputs:
    %   battery_GD: Struct containing the current battery data.
    %
    % Outputs:
    %   battery_GD: Updated battery structure after automatic design.
    %   SysSpec: Updated system specifications.

    % Extract SysSpec for convenience
    SysSpec = battery_GD.SysSpec;

    % Step 1: Create a GUI for user inputs
    inputFig = figure('Name', 'Automatic Design Inputs', 'NumberTitle', 'off', ...
                      'Position', [500, 300, 400, 300], 'MenuBar', 'none', 'ToolBar', 'none');

    % Input fields with labels
    uicontrol(inputFig, 'Style', 'text', 'Position', [50, 250, 150, 25], ...
              'String', 's_mod_max:', 'HorizontalAlignment', 'right');
    s_mod_maxInput = uicontrol(inputFig, 'Style', 'edit', 'Position', [210, 250, 120, 25], ...
                               'String', num2str(SysSpec.s_mod_max), ...
                               'Callback', @(src, ~) updateSysSpec('s_mod_max', src.String));

    uicontrol(inputFig, 'Style', 'text', 'Position', [50, 200, 150, 25], ...
              'String', 'U_mod_nom (V):', 'HorizontalAlignment', 'right');
    U_mod_nomInput = uicontrol(inputFig, 'Style', 'edit', 'Position', [210, 200, 120, 25], ...
                               'String', num2str(SysSpec.U_mod_nom), ...
                               'Callback', @(src, ~) updateSysSpec('U_mod_nom', src.String));

    uicontrol(inputFig, 'Style', 'text', 'Position', [50, 150, 150, 25], ...
              'String', 'I_mod_max (A):', 'HorizontalAlignment', 'right');
    I_mod_maxInput = uicontrol(inputFig, 'Style', 'edit', 'Position', [210, 150, 120, 25], ...
                               'String', num2str(SysSpec.I_mod_max), ...
                               'Callback', @(src, ~) updateSysSpec('I_mod_max', src.String));

    uicontrol(inputFig, 'Style', 'text', 'Position', [50, 100, 150, 25], ...
              'String', 'C_mod_min (Ah):', 'HorizontalAlignment', 'right');
    C_mod_minInput = uicontrol(inputFig, 'Style', 'edit', 'Position', [210, 100, 120, 25], ...
                               'String', num2str(SysSpec.C_mod_min), ...
                               'Callback', @(src, ~) updateSysSpec('C_mod_min', src.String));

    % Proceed button to finalize inputs
    uicontrol(inputFig, 'Style', 'pushbutton', 'Position', [150, 30, 100, 40], ...
              'String', 'Save', 'Callback', @saveCallback);

    % Block execution until GUI is closed
    uiwait(inputFig);

    % Nested functions for updating SysSpec and closing the GUI
    function updateSysSpec(field, value)
        SysSpec.(field) = str2double(value);
        battery_GD.SysSpec = SysSpec;
    end

    function saveCallback(~, ~)
        % Resume execution and close the GUI
        uiresume(inputFig);
        delete(inputFig);
    end
end