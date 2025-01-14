%% Single Track Model
% Input driving cycle:
%   time    [s]
%   speed   [m/s]
%   gear    []
% Output tire:
%   Torque  [Nm]
%   RPM     [1/min]
%% Author
% Frederic Brenner

driving_cycle = "WLTP_class_3";
load(".\data\DrivingCycles\"+driving_cycle+".mat");

[torque, rpm] = calculate_tire_speed(dc.time, dc.speed, config_vehicle)

% plot the driving cycle
plot(rpm, abs(torque),'.')
ylabel("T_{tire} [Nm]")
xlabel("RPM_{tire} [min^{-1}]")
% ylim([-400 400])
% xlim([0 15000])
hold on
