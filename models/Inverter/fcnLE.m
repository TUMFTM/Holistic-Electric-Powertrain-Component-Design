udc = config_battery.U_sys_nom;                     %battery DC Voltage                    (V)
fswitch = 10000;                                    %switch frequency                      (Hz)
%% IGBT: Equivalent circuit thermal model RC
R1igbt = 0.0062;
tau1igbt = 0.0005;
C1igbt = tau1igbt/R1igbt;
R2igbt = 0.022;
tau2igbt = 0.02;
C2igbt = tau2igbt/R2igbt;
R3igbt = 0.0238;
tau3igbt = 0.058;
C3igbt = tau3igbt/R3igbt;
R4igbt = 0.038;
tau4igbt = 0.45;
C4igbt = tau4igbt/R4igbt;
R5igbt = 0.01;
tau5igbt = 2.19;
C5igbt = tau5igbt/R5igbt;
R6igbt = 0.1;
%
R1diode = 0.0126;
tau1diode = 0.0006;
C1diode = tau1diode/R1diode;
R2diode = 0.0391;
tau2diode = 0.017;
C2diode = tau2diode/R2diode;
R3diode = 0.0417;
tau3diode = 0.053;
C3diode = tau3diode/R3diode;
R4diode = 0.0331;
tau4diode = 0.39;
C4diode = tau4diode/R4diode;
R5diode = 0.0135;
tau5diode = 3.2;
C5diode = tau5diode/R5diode;
%%
EV_array = zeros(1,1);              %consumption of BEV in different driving cycles
Efficiency_array = zeros(1,1);      %efficiency of inverter in different driving cycles
Conduction_loss_array = zeros(1,1); %conduction loss of inverter in different driving cycles
Switching_loss_array = zeros(1,1);  %switching loss of inverter in different driving cycles
Inverter_loss_array = zeros(1,1);   %total loss of inverter in different driving cycles
% All energy values are converted to kWh/100km
k = 0;
efficiency_calc = 'Main_EfficiencyCalculation_Motor_Revised.m';
%%
% model = 'Main_Motor_Revised_SiC.slx';              % for SiC MOSFET inverter
model = 'Main_ComplexTherma_Test_IGBT_Revised.slx';  % for IGBT inverter
%%
driving_cycle = "WLTP_class_3_extra_high";
load("data/DrivingCycles/"+driving_cycle+".mat");
DRIVING_CYCLE = [dc.time,dc.speed];
k = k + 1;
time = length(DRIVING_CYCLE);
%%
sim(model);
run(efficiency_calc);
EV_array(k) = A+S+M;
Efficiency_array(k) = 1-S/(S+B);
Conduction_loss_array(k) = Conduction_loss;
Switching_loss_array(k) = Switching_loss;
Inverter_loss_array(k) = S;