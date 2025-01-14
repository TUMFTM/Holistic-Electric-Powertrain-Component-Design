function [xC, yC, zC] = calc_motor_isolines(mat_speed,mat_torque,mat_efficiency)
%CALC_MOTOR_ISOLINES Summary of this function goes here
%   Detailed explanation goes here

    % mat_speed = config_motor.Speed;
    % mat_torque = config_motor.Shaft_Torque;
    % mat_efficiency = config_motor.Efficiency;

    figure('Name', 'Motor efficiency isolines');
    K = contourf(mat_speed, mat_torque, mat_efficiency, levels=10);
    
    iC = 1;
    idxC = 1;
    hold on
    while iC < size(K,2)
      nP = K(2,iC); % number of points in current contour
      zC{idxC} = K(1,iC);
      xC{idxC} = K(1,iC+(1:nP)); % x coordinates of current contour
      yC{idxC} = K(2,iC+(1:nP)); % y coordinates of current contour
      iC = iC+nP+1;    % Start-point of next contour
      idxC = idxC + 1; % next contourline index
      plot(xC{idxC-1},yC{idxC-1},'b.-') %% If you want to look at the separate
    end

    close
end
