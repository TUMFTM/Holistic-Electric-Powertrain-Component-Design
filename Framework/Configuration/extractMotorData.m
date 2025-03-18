function [frontMotorData, rearMotorData] = extractMotorData(dataTable, type)
    % Determine the search term based on the motor type
    switch type
        case 'PSM'
            frontMotorType = 'Permanent Magnet / AC Synchronous';
            rearMotorType = 'Permanent Magnet/AC Synchronous';
        case 'ASM'
            frontMotorType = 'Induction / AC Asynchronous';
            rearMotorType = 'Induction / AC Asynchronous';
        otherwise
            error('Invalid motor type. Please specify either ''PSM'' or ''ASM''.');
    end

    % Initialize conditions for front and rear motors
    frontCondition = false(height(dataTable), 1); % Preallocate logical array for front condition
    
    % Loop through each row to check the conditions for FrontE_Motor1Type
    for i = 1:height(dataTable)
        if ~isempty(dataTable.FrontE_Motor1Type{i}) && isequal(dataTable.FrontE_Motor1Type{i}, frontMotorType)
            frontCondition(i) = true;
        end
    end

    % Initialize conditions for rear motors
    rearCondition = false(height(dataTable), 1);
    for i = 1:height(dataTable)
        currentValue = dataTable.RearE_Motor1Type{i};
        if ~isempty(currentValue) && isequal(currentValue, rearMotorType)
            rearCondition(i) = true;
        end
    end

    % Extract rows and columns for front motor
    frontMotorData = dataTable(frontCondition, {'FrontE_Motor1Power_kW_', ...
                                                 'FrontE_Motor1Torque_Nm_', ...
                                                 'TotalWeight_kg__FrontE_motor_'});
                                             
    % Extract rows and columns for rear motor
    rearMotorData = dataTable(rearCondition, {'RearE_Motor1Power_kW_', ...
                                              'RearE_Motor1Torque_Nm_', ...
                                              'TotalWeight_kg__RearE_motor_'});
end