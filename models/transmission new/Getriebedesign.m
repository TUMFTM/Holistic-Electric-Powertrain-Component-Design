function [Gearbox] = Getriebedesign(M,iG,AnzahlStufen)
% Function for the design of a parallel-axis gearbox, single- or two-stage
% Input: allowable torque M,
%        total gear ratio iG,
%        number of stages AnzahlStufen (1 or 2)

Mmax = M;
iges = iG;

%% Reference values from A2Mac1 for VW ID.3 and Tesla Model 3
% ID.3 gearbox from A2Mac1 (i1 = 2.95)  (i2 = 3.9)
% Stage 1: Pinion: z=23   dk=50mm   b=35mm
%          Gear: z=68   dk=135mm   b=32mm
%
% Stage 2: Pinion: z=20   dk=63mm   b=49mm
%          Gear: z=78   dk=217mm   b=45mm

% Tesla Model 3 gearbox from A2Mac1 (i1 = 2.61)  (i2 = 3.46)
% Stage 1: Pinion: z=31   dk=58mm   b=37mm
%          Gear: z=81   dk=142mm   b=37mm
%
% Stage 2: Pinion: z=24   dk=67mm   b=46mm
%          Gear: z=83   dk=212mm   b=46mm


%% Calculations

simga_H_lim = 1500; %N/mm^2, Endurance/Time strength value from Naunheimer, Fahrzeuggetriebe, 2019, Tab. 7.1

alpha = 20; % Degrees; standard pressure angle; 20° = standard value
beta = [23.4 18.4]; % Degrees; helix angle; adopted from Gao

% Literature recommendations for the number of stages
% Single-stage up to i=5, otherwise two-stage according to Sauer, Konstruktionselemente, 2018, p.480
% Single-stage up to i=6, otherwise two-stage according to Niemann, Getriebe allgemein, 2003, p.261

if AnzahlStufen > 1 % Two-stage, requires distribution of the gear ratio between stage 1 and stage 2
    
    % Distribution of the gear ratio between stage 1 and stage 2 according to Niemann, Getriebe allgemein, 2003, p.261
    % i(1) = 0.8*iges^(2/3);
    % i(2) = iges/i(1);
    
    % Distribution of the gear ratio between stage 1 and stage 2 according to Sauer,
    % Konstruktionselemente des Maschinenbaus 2, p.480 for 5<iges<15
    % based on Römhild Iris 1993
    i(1) = 0.7332*iges^0.6438;
    i(2) = iges/i(1);
    
    % Considering the integrated differential in the last stage may require
    % a larger gear ratio in stage 2 to create space for the differential
    % The following lines can be activated as an overwrite:
    % i(2) = [3.9]; %Overwrite ID3 [2.95 3.9]
    % i(2) = [3.46]; %Overwrite Tesla M3 [2.61 3.46]
    i(1) = iges/i(2);

    % Empirical values for pinion teeth number from Niemann, Getriebe allgemein, 2003, p.269, case-hardened
    i_inter = [1 2 4 8];
    z_inter1 = [32 29 25 22]; % high speed
    z_inter2 = [26.5 24 20.5 18]; % medium speed
    
    % Interpolation of teeth number for stage 1 and 2 based on empirical values
    z(1) = round(interp1(i_inter,z_inter1,i(1),'linear'),0); % high speeds
    z(2) = round(interp1(i_inter,z_inter2,i(2),'linear'),0); % medium speeds
    
    % To overwrite with known teeth numbers, activate:
    % z = [23 20]; % ID.3 pinion teeth number;
    % z = [31 24]; % Tesla Model 3 input

else % Single-stage case
    i = iges;
    if iges > 6 % Warning if the recommended gear ratio range is exceeded
        disp("Caution, total ratio too big for one stage according to literature!")
    end
    
    % Empirical values for pinion teeth number from Niemann, Getriebe allgemein, 2003, p.269, case-hardened
    i_inter = [1 2 4 8];
    z_inter1 = [32 29 25 22]; % high speed

    % Interpolation of pinion teeth number based on empirical values
    z(1) = round(interp1(i_inter,z_inter1,i(1),'linear','extrap'),0); % high speeds
    
    % In the single-stage case, only one helix angle is necessary
    beta = beta(1);
end

% % Calculation of gear teeth number(s)
z2 = round(z.*i);


%% Shaft diameter
% Shaft diameter required for seal size and bearing selection

% Equivalent moment due to bending and torsion loading: Factor 1.2 ... 2.5
% Source: Niemann Winter Maschinenelemente 1, 2019, Formula 17.10
% Mv = 1.5*Mmax; -> used below

% Steel assumption: 16MnCr5 sigma_bW = 520 N/mm^2 from Naunheimer, Fahrzeuggetriebe, 2019, p.478
% Safety factor 1.5, size factor 0.85, surface factor 0.8

sigma_b_zul = 520*0.85*0.8/1.5; % allowable bending stress according to Naunheimer, Fahrzeuggetriebe, 2019, p.478

for j=1:AnzahlStufen+1
    Mv = 1.5*Mmax*prod(i(1:j-1)); % Equivalent moment
    if j > AnzahlStufen % if output shaft
        Mv = Mv/2; % Split the moment between left and right
    end
    d(j) = 2.17 * nthroot(Mv*1e3/sigma_b_zul,3); % Shaft diameter according to Naunheimer, Fahrzeuggetriebe, 2019, p.476
end


%% Center distance a, tooth width b, and normal module mn

M(1) = Mmax; % Torque stage 1
M(2) = M(1)*i(1); % Torque stage 2
for j = 1:AnzahlStufen
    % Center distance according to Naunheimer, Fahrzeuggetriebe, 2019, p.288 
    a(j) = 0.255 * nthroot((M(j)*1e3*(i(j)+1)^4)/i(j),3);
    
    % Pitch circle diameter according to Niemann, Getriebe allgemein, 2003, p.276
    d1(j) = 2*a(j)/(1+z2(j)/z(j));
    d2(j) = 2*a(j)/(1+z(j)/z2(j));
    
    % Tooth width according to Naunheimer, Fahrzeuggetriebe, 2019, p.289
    b(j) = 4278e2 * 0.65 * (M(j)*1e3*(i(j)+1))/(d1(j)^2*i(j)*simga_H_lim^2);

    % Normal module according to Wittel, Roloff/Matek Maschinenelemente, 2021, p.813
    mn(j) = (2*a(j)*cosd(beta(j)))/(z(j)+z2(j));

    disp(['Pitch circle diameter stage ',int2str(j),': ',sprintf('%0.2f  %0.2f',d1(j),d2(j)),' mm'])
    disp(['Tooth width stage ',int2str(j),': ',sprintf('%0.2f',b(j)),' mm'])
end


%% Bearing selection

% SKF bearing catalog series 160: 16002 to 16020
di = [15 17 20:5:100]; % Inner diameter
da = [32 35 42 47 55 62 68 75 80 90 95 100 110 115 125 130 140 145 150]; % Outer diameter
c_stat = [2.85 3.25 4.05 4.75 7.35 8.15 10.2 10.8 11.4 14 15 19.6 25 27 31.5 33.5 39 41.5 44]*1e3; % Static load rating

for j = 1:AnzahlStufen+1 % Iterate over shafts, both bearings per shaft identical
    % Select the next larger bearing
    Lager_index = sum(d(j) > di)+1;
    % Assign values from bearing catalog
    dLager(j) = di(Lager_index);
    DLager(j) = da(Lager_index);
    C0Lager(j) = c_stat(Lager_index);
end


%% Saving the design in the Gearbox struct
Gearbox.Stufen = AnzahlStufen; % Number of stages
Gearbox.iges = iges; % Total gear ratio
Gearbox.alpha = alpha; % Standard pressure angle
Gearbox.beta = beta; % Helix angle
Gearbox.z1 = z; % Number of teeth on pinion
Gearbox.z2 = z2; % Number of teeth on gear(s)
Gearbox.d1 = d1; % Pitch circle diameter of pinion
Gearbox.d2 = d2; % Pitch circle diameter of gear(s)
Gearbox.a = a; % Center distance(s)
Gearbox.b = b; % Tooth width(s)
Gearbox.mn = mn; % Normal module(s)
Gearbox.d = d; % Shaft diameter
Gearbox.dLager = dLager; % Bearing inner diameter
Gearbox.DLager = DLager; % Bearing outer diameter
Gearbox.C0Lager = C0Lager; % Static load ratings of bearings

end

