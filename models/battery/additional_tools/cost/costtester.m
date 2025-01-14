%% costtester v0.1
% This interface for sim_BTMS was created as part of Global Drive 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                           Versions                                    %
% v0.1 - first setup                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        known Problems                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init
clear
clc

%% Set Testparameters

h.n = 2;                % Number of runs 
h.E_sys_max = 100;      % Maximal Energy in kWh

% Possible inputs
p.cellchemistry = ["NMC_111";"NMC_532";"NMC_622";"NMC_811";"NCA";"LFP";"LMO";"Natrium"];
p.E_sys = [10 h.E_sys_max]; %[kWh]

t.E_sys = round(rand(h.n,1) * (p.E_sys(2)-p.E_sys(1)) + p.E_sys(1)); %[kWh]
t.cellchemistry = (round(rand(h.n,1).*numel(p.cellchemistry)));

%% Test
for i = 1: h.n
    % try 
        % Cost calculated in Excel using ActivX
        A = costdata_evaluator_testversion(t.E_sys(i),p.cellchemistry(t.cellchemistry(i)));

        configs_6_BTMS_passed.BatPara.cellchemistry = p.cellchemistry(t.cellchemistry(i));
        configs_6_BTMS_passed.SysInfo.E_sys = t.E_sys(i);
        B = costcode_bat(configs_6_BTMS_passed);
    % catch
    % end
end
      

