% Preprocess motormap to extract maximum torque for each speed
motormap = evalin("base","motormap");

secondMotor = evalin("base","secondMotor");

max_torque_map = max(motormap.Shaft_Torque, [], 2); % Maximum torque per row
speed_vector = motormap.Speed(:, 1); %  Assuming all rows in motormap.Speed are identical
efficiency_vector = motormap.Efficiency(:, 1)/100; % Extract the first column of efficiency
lossData = motormap.Total_Loss(:, 1); 

% Define torque and speed vectors (assumed to be provided)
%max_torque_map = max(motormap2.Shaft_Torque, [], 2); % Max torque per row
%speed_vector = motormap2.Speed(:, 1); % Assuming all rows have identical speed values


if secondMotor
    motormap2 = evalin("base","motormap2");
    max_torque_map2 = max(motormap2.Shaft_Torque, [], 2); % Maximum torque per row
    speed_vector2 = motormap2.Speed(:, 1); % Assuming all rows in motormap2.Speed are identical
    efficiency_vector2 = motormap2.Efficiency(:, 1)/100; % Extract the first column of efficiency

    % Define new speed and torque vectors for interpolation
    %speedVec2 = linspace(min(speed_vector2), max(speed_vector2), 150); % Generate 150 points for speed
    %torqueVec2 = linspace(min(max_torque_map2), max(max_torque_map2), 120); % Generate 120 points for torque
    lossData2 = motormap2.Total_Loss(:, 1); % Maximum loss per row
end