function [battery_GD] = sustainability_estimator(battery_GD)

%% sustainability_estimator v0.1
% This Code was created as part of Global Drive 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                           Versions                                    %
% v0.1 - first setup of System Parametrisation                          %
                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        known Problems                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load database
ReCiPe = load("ReCiPe_2016_1kWh.mat");
ReCiPe = ReCiPe.ReCiPe_2016_1kWh;

if any(battery_GD.BatPara.cellchemistry == ReCiPe.cellchemistry)
    h.pos = find(ismember(ReCiPe.cellchemistry, battery_GD.BatPara.cellchemistry));
    for i = 1:numel(h.pos)
        sust.(string(ReCiPe.country(h.pos(i)))) = ReCiPe(h.pos(i),5:size(ReCiPe, 2)).*battery_GD.SysInfo.E_sys;
        sust.calculationmethode = "ReCiPe_2016";
    end

else
    warning("No Sustainability data for this cellchemistry available")
    sust.not_available = 1;
end



battery_GD.SysInfo.Sustainability = sust;
end