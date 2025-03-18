function  [f, fields] = create_parameter_interaction_panel()

       % Create main figure
    %f = uifigure('Name', 'Vehicle Design Parameters', 'Position', [100, 100, 1000, 800]);
    
    screenSize = get(0, 'ScreenSize'); % Get screen size [left, bottom, width, height]
    figWidth = 1400;  % Set desired width
    figHeight = 900;  % Set desired height
    figX = (screenSize(3) - figWidth) / 2;  % Center X position
    figY = (screenSize(4) - figHeight) / 2; % Center Y position
    
    f = uifigure('Name', 'Vehicle Design Parameters', ...
                 'Position', [figX, figY, figWidth, figHeight]);

    % Initialize panel dimensions and spacing
    panelWidth = 400;
    panelHeight = 380;
    panelSpacing = 50;
    
    fieldWidth = 200;
    fieldHeight = 22;
    fieldSpacingY = 30;
    labelWidth = 150;
    labelHeight = 22;
    startX = 10;
    startY = panelHeight - 100;
    
    
    % ---- Battery Parameters Panel ----
    batteryPanel = uipanel(f, 'Title', 'Battery Parameters', ...
        'Position', [50, 520, panelWidth, 380]);
    currentY = startY;
    
    % Battery Parameters
    batteryChemistryLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Battery Chemistry');
    batteryChemistryField = uidropdown(batteryPanel, 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], ...
        'Items', {'Cylindrical NMC811 5Ah', 'Prismatic NMC811 168Ah',  'Cylindrical NMC721 27Ah',  'Pouch NMC811 109Ah',  'Prismatic LFP 249Ah', 'Prismatic NMC811 529Ah'}, 'Value', 'Cylindrical NMC811 5Ah');
    currentY = currentY - fieldSpacingY;
    
    serialConnectionsLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Serial Connections');
    serialConnectionsField = uieditfield(batteryPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    parallelConnectionsLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Parallel Connections');
    parallelConnectionsField = uieditfield(batteryPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    % Assign variables for each battery parameter
    batteryEnergyLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Battery Energy (Wh)');
    batteryEnergyField = uieditfield(batteryPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    batteryVoltageLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Battery Voltage (V)');
    batteryVoltageField = uieditfield(batteryPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    batteryCapacityLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Battery Capacity (Ah)');
    batteryCapacityField = uieditfield(batteryPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    batteryCurrentLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Battery Current (A)');
    batteryCurrentField = uieditfield(batteryPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    batteryPowerLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Battery Power (W)', 'Visible', 'off');
    batteryPowerField = uieditfield(batteryPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0, 'Visible', 'off');

    batteryMassLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Battery Mass (kg)');
    batteryMassField = uieditfield(batteryPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    batteryVolumeLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Battery Volume (L)');
    batteryVolumeField = uieditfield(batteryPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    
    currentY = currentY - fieldSpacingY;

    batteryEfficiencyLabel = uilabel(batteryPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Battery Efficiency (0-1)');
    batteryEfficiencyField = uieditfield(batteryPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    
    % ---- Cell Parameters Panel ----
    cellPanel = uipanel(f, 'Title', 'Cell Parameters', ...
        'Position', [50, 225, panelWidth, 250]);
    currentY = 200;
    
    cellEnergyLabel = uilabel(cellPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Cell Energy (Wh)');
    cellEnergyField = uieditfield(cellPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    cellVoltageLabel = uilabel(cellPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Cell Voltage (V)');
    cellVoltageField = uieditfield(cellPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    cellCapacityLabel = uilabel(cellPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Cell Capacity (Ah)');
    cellCapacityField = uieditfield(cellPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    cellCurrentLabel = uilabel(cellPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Cell Current (A)');
    cellCurrentField = uieditfield(cellPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    cellMassLabel = uilabel(cellPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Cell Mass (kg)');
    cellMassField = uieditfield(cellPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    cellVolumeLabel = uilabel(cellPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Cell Volume (m^3)');
    cellVolumeField = uieditfield(cellPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    
    
    % ---- Motor Parameters Panel ----
    motorPanel = uipanel(f, 'Title', 'Motor Parameters', ...
        'Position', [500, 500, panelWidth, panelHeight]);
    currentY = startY;
    
    motorPeakPowerLabel = uilabel(motorPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Motor Peak Power (W)');
    motorPeakPowerField = uieditfield(motorPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    motorPeakTorqueLabel = uilabel(motorPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Motor Peak Torque (Nm)');
    motorPeakTorqueField = uieditfield(motorPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;

    motorWeightLabel = uilabel(motorPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Motor Weight (kg)');
    motorWeightField = uieditfield(motorPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    motorEfficiencyLabel = uilabel(motorPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Motor Efficiency (0-1)');
    motorEfficiencyField = uieditfield(motorPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    
    % Add Motor Type Label
    motorTypeLabel = uilabel(motorPanel, ...
        'Position', [startX, currentY - 50, labelWidth, labelHeight], ...
        'Text', 'Motor Type');
    
    % Add Motor Type Dropdown
    motorTypeField = uidropdown(motorPanel, ...
        'Position', [startX + labelWidth + 10, currentY - 50, fieldWidth, fieldHeight], ...
        'Items', {'PSM', 'ASM'}, ...
        'Value', 'PSM', ... % Default value
        'Editable', 'off', ... % Dropdown is not editable
        'Tooltip', 'Select the motor type (PSM or ASM)');
    
    
    % ---- Gearbox Parameters Panel ----
    gearboxPanel = uipanel(f, 'Title', 'Gearbox Parameters', ...
        'Position', [50, -50, panelWidth, 250]);
    currentY = 250 - (1 * (fieldHeight + fieldSpacingY));
    %currentY = startY;
    
    gearboxRatioLabel = uilabel(gearboxPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Gearbox Ratio');
    gearboxRatioField = uieditfield(gearboxPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    gearboxMaxTorqueLabel = uilabel(gearboxPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Gearbox Max Torque (Nm)');
    gearboxMaxTorqueField = uieditfield(gearboxPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    gearboxWeightLabel = uilabel(gearboxPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Gearbox Weight (kg)');
    gearboxWeightField = uieditfield(gearboxPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    currentY = currentY - fieldSpacingY;
    
    gearboxEfficiencyLabel = uilabel(gearboxPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Gearbox Efficiency (0-1)');
    gearboxEfficiencyField = uieditfield(gearboxPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    
    
    % ---- Inverter Parameters Panel ----
    inverterPanel = uipanel(f, 'Title', 'Inverter Parameters', ...
        'Position', [500, 100, panelWidth, panelHeight]);
    currentY = startY;
    
    inverterEfficiencyLabel = uilabel(inverterPanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Inverter Efficiency');
    inverterEfficiencyField = uieditfield(inverterPanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0);
    
    % Add Inverter Type Label
    inverterTypeLabel = uilabel(inverterPanel, ...
        'Position', [startX, currentY - 50, labelWidth, labelHeight], ...
        'Text', 'Inverter Type');
    
    % Add Inverter Type Dropdown Field
    inverterTypeField = uidropdown(inverterPanel, ...
        'Position', [startX + labelWidth + 10, currentY - 50, fieldWidth, fieldHeight], ...
        'Items', {'MOSFET', 'IGBT'}, ... % Dropdown options
        'Value', 'MOSFET', ... % Default value
        'Tooltip', 'Select the inverter type (MOSFET or IGBT)');
    % ---- Vehicle Performance Parameters Panel ----
    vehiclePanel = uipanel(f, 'Title', 'Vehicle Parameters', ...
        'Position', [950, 100, panelWidth, panelHeight]);
    currentY = startY;
    
    vehicleMaxSpeedLabel = uilabel(vehiclePanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Vehicle Max Speed (km/h)');
    vehicleMaxSpeedField = uieditfield(vehiclePanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0, 'Editable', 'off');
    currentY = currentY - fieldSpacingY;
    
    vehicleAccTimeLabel = uilabel(vehiclePanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Acceleration Time (s)');
    vehicleAccTimeField = uieditfield(vehiclePanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0, 'Editable', 'off');
    currentY = currentY - fieldSpacingY;
    
    vehicleConsumptionLabel = uilabel(vehiclePanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Vehicle Consumption (kWh/100km)');
    vehicleConsumptionField = uieditfield(vehiclePanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0,'Editable', 'off' );
    currentY = currentY - fieldSpacingY;
    
    vehicleRangeLabel = uilabel(vehiclePanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Vehicle Range (km)');
    vehicleRangeField = uieditfield(vehiclePanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0, 'Editable', 'off');
    currentY = currentY - fieldSpacingY;
    
    vehicleEfficiencyLabel = uilabel(vehiclePanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Vehicle Efficiency');
    vehicleEfficiencyField = uieditfield(vehiclePanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0, 'Editable', 'off');
    currentY = currentY - fieldSpacingY;
    
    vehicleMassLabel = uilabel(vehiclePanel, 'Position', [startX, currentY, labelWidth, labelHeight], 'Text', 'Vehicle Mass (kg)');
    vehicleMassField = uieditfield(vehiclePanel, 'numeric', 'Position', [startX + labelWidth + 10, currentY, fieldWidth, fieldHeight], 'Value', 0, 'Editable', 'off');
    %outputFields = struct(motorMassField);
    
    
    % Map parameter names to their UI components (assuming UI components are defined in the main function)
        fields = struct(...
        'battery_energy', batteryEnergyField, ...
        'battery_voltage', batteryVoltageField, ...
        'battery_capacity', batteryCapacityField, ...
        'battery_current', batteryCurrentField, ...
        'battery_mass', batteryMassField, ...
        'battery_volume', batteryVolumeField, ...
        'battery_power', batteryPowerField, ...
        'battery_efficiency', batteryEfficiencyField, ...
        'cell_energy', cellEnergyField, ...
        'cell_voltage', cellVoltageField, ...
        'cell_capacity', cellCapacityField, ...
        'cell_current', cellCurrentField, ...
        'cell_mass', cellMassField, ...
        'cell_volume', cellVolumeField, ...
        'serial_conections', serialConnectionsField, ...
        'parallel_connections', parallelConnectionsField, ...
        'motor_peak_power', motorPeakPowerField, ...
        'motor_peak_torque', motorPeakTorqueField, ...
        'motor_mass', motorWeightField, ...
        'motor_efficiency', motorEfficiencyField, ...
        'gearbox_ratio', gearboxRatioField, ...
        'gearbox_max_torque', gearboxMaxTorqueField, ...
        'gearbox_mass', gearboxWeightField, ...
        'gearbox_efficiency', gearboxEfficiencyField, ...
        'inverter_efficiency', inverterEfficiencyField, ...
        'vehicle_max_speed', vehicleMaxSpeedField, ...
        'vehicle_acc_time', vehicleAccTimeField, ...
        'vehicle_consumption', vehicleConsumptionField, ...
        'vehicle_range', vehicleRangeField, ...
        'vehicle_efficiency', vehicleEfficiencyField, ...
        'vehicle_mass', vehicleMassField, ...
        'battery_chemistry', batteryChemistryField, ...
        'motor_type', motorTypeField,...
        'inverter_type', inverterTypeField);
     
        % Set non-editable and visually distinct properties for the specified fields
    nonEditableFields = { ...
        fields.cell_energy, ...
        fields.cell_voltage, ...
        fields.cell_capacity, ...
        fields.cell_current, ...
        fields.cell_mass, ...
        fields.cell_volume, ...
        fields.motor_efficiency, ...
        fields.vehicle_max_speed, ...
        fields.vehicle_acc_time, ...
        fields.vehicle_consumption, ...
        fields.vehicle_range, ...
        fields.vehicle_efficiency, ...
        fields.vehicle_mass, ...
    };
    
    % Apply non-editable and styling properties to all fields in the list
    for i = 1:length(nonEditableFields)
        field = nonEditableFields{i};
        %field.Editable = 'off';                % Make non-editable
        field.BackgroundColor = [0.9, 0.9, 0.9]; % Light gray background
        field.FontWeight = 'bold';            % Bold text
    end


end