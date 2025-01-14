function [battery] = configs_elector(battery,configs_6_BTMS_passed)

%% configs_elector v0.1
% This function is part of the Global Drive 2023 project
% It chooses best fitting battery configuration of the ones found by
% sim_BTMS_1_system_setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                           Versions                                    %
% v0.1 - creation                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        known Problems                                 %
% - none                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

configs_values = Inf (size(configs_6_BTMS_passed,2),1); % only mass so far
for i = 1:size(configs_6_BTMS_passed,2)
    configs_values(i) = configs_6_BTMS_passed(i).SysInfo.mass_sys;
end

% Add weighting of different cellparameters here

[~,min_idx]=min(configs_values);

battery = configs_6_BTMS_passed(min_idx);

fprintf("Successfully chose one battery configuration\n")

end
