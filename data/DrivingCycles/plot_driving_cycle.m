% load('WLTP_class_3.mat');
% load('USA_FTP_75.mat');
load('EUROPE_ARTEMIS_MOTORWAY.mat');

figure('units','centimeters','position',[0,0,18,10])
set(gca,'fontname','Arial') 
set(gca,'FontSize',11)

plot(dc.time, dc.speed.*3.6, 'k', 'LineWidth', 1)

xlabel('Time in s', 'FontName', 'Arial')
ylabel('Velocity in km/h', 'FontName', 'Arial')

if false
    lines = [0, 589, 589+432, 589+432+454];
    xl = xline(lines, 'b-.', {'Low', 'Medium', 'High', 'Extra High'}, 'FontName', 'Arial');
    % xl(1).LabelHorizontalAlignment = 'left';
    % xl(2).LabelHorizontalAlignment = 'left';
    % xl(3).LabelHorizontalAlignment = 'left';
    % xl(4).LabelHorizontalAlignment = 'left';
    xl(1).LabelOrientation = 'horizontal';
    xl(2).LabelOrientation = 'horizontal';
    xl(3).LabelOrientation = 'horizontal';
    xl(4).LabelOrientation = 'horizontal';
end

