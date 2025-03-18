function [motorInfo,  img] = displayMotorInformation(motormapName, config_motor, second)
    % DISPLAYMOTORINFORMATION - Extracts and displays motor specifications 
    % and associated images for the selected motor.
    %
    % Inputs:
    %   motormapName - Name of the motor map file (string)
    %
    % Outputs:
    %   motorInfo - Table containing extracted motor parameters
    %   img       - Image of the motor (if available)

    
    % Check if config_motor structure exists in the workspace
    if exist('config_motor', 'var') && isstruct(config_motor)
        if second
        motorInfo1 = table( ...
            {'Motor Type'; 'Voltage (V)'; 'Number of Motors'; 'Peak Torque (Nm)'; 'Peak Power (kW)'; ...
             'Motor Current (Arms)'; 'Total Peak Torque (Nm)'; 'Total Peak Power (kW)'}, ...
            {config_motor.type2; config_motor.voltage2; config_motor.motorsCount2; ...
             config_motor.torque2; config_motor.power2; config_motor.current2; ...
             config_motor.totalTorque2; config_motor.totalPower2}, ...
            'VariableNames', {'Parameter', 'Value'} ...
        );
        else    
        % Define motor information fields for Motor 1
        motorInfo1 = table( ...
            {'Motor Type'; 'Voltage (V)'; 'Number of Motors'; 'Peak Torque (Nm)'; 'Peak Power (kW)'; ...
             'Motor Current (Arms)'; 'Total Peak Torque (Nm)'; 'Total Peak Power (kW)'}, ...
            {config_motor.type1; config_motor.voltage1; config_motor.motorsCount1; ...
             config_motor.torque1; config_motor.power1; config_motor.current1; ...
             config_motor.totalTorque1; config_motor.totalPower1}, ...
            'VariableNames', {'Parameter', 'Value'} ...
        );
        end
        
        % Initialize motorInfo output
        motorInfo = motorInfo1;
        
        

    else
        % Display a warning if motor configuration is missing
        disp('Motor configuration data (config_motor) is missing or invalid.');
        motorInfo = table(); % Return empty table if no valid data
    end

    % Attempt to load and display motor image
    img = [];
    try
        imageFile = replace(motormapName, '.mat', '.png'); % Convert .mat to .png
        img = imread(imageFile);
        % imshow(img); % Uncomment if running in a live script environment
    catch
        warning('Image file could not be found or loaded.');
    end
end