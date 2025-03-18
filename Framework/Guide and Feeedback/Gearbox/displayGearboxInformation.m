function gearbox_info_table = displayGearboxInformation(GB, GB_output)
    % Initialize data arrays for the table
    parameter = {};
    value = {};
    
    % Add stage-specific information
    for j = 1:GB.Stufen
        parameter{end+1} = ['Pitch Circle Diameter Stage ', int2str(j)];
        value{end+1} = sprintf('%0.2f  %0.2f mm', GB.d1(j), GB.d2(j));
        
        parameter{end+1} = ['Tooth Width Stage ', int2str(j)];
        value{end+1} = sprintf('%0.2f mm', GB.b(j));
    end
    
    % Add general gearbox information
    parameter{end+1} = 'Gearbox Weight (Kg)';
    value{end+1} = sprintf('%0.2f', GB.m_gearbox);
    
    parameter{end+1} = 'Gearbox Volume (m^3)';
    value{end+1} = sprintf('%0.2f', GB.V_gearbox);
    
    parameter{end+1} = 'Gear Losses Power (W)';
    value{end+1} = sprintf('%0.2f', GB_output.gear_loss);
    
    parameter{end+1} = 'Seal Losses Power (W)';
    value{end+1} = sprintf('%0.2f', GB_output.seal_loss);
    
    parameter{end+1} = 'Bearing Losses Power (W)';
    value{end+1} = sprintf('%0.2f', GB_output.bearing_loss);
    
    parameter{end+1} = 'Total Losses Power (W)';
    value{end+1} = sprintf('%0.2f', GB_output.total_loss);

    % Create a table
    gearbox_info_table = table(parameter', value', ...
        'VariableNames', {'Parameter', 'Value'});
    
    % Display the table inline in the Live Script
    %disp('Gearbox Information:');
    %disp(gearbox_info_table);
end