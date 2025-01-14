function [Gearbox] = Getriebedesign_koaxial(M,iG)
% Function for the design of a coaxial gearbox (two-stage)
% Input: allowable torque M, 
%        total gear ratio iG, 
Mmax = M;
iges = iG;


%% Calculation

simga_H_lim = 1500; % N/mm^2, endurance/time strength value from Naunheimer, Fahrzeuggetriebe, 2019, Table 7.1

alpha = 20; % Degrees; standard pressure angle; 20° = standard value
beta = [23.4 18.4]; % Degrees; helix angle; adopted from Gao
AnzahlStufen = 2; % Must be two-stage
    
% Distribution of the gear ratio across both stages, such that a(1) == a(2) for coaxiality
% Formula for center distance according to Naunheimer, Fahrzeuggetriebe, 2019, p.288
syms i1; % symbolic variable for gear ratio stage 1
assume(i1,'real');
assumeAlso(i1>0);
M(1) = Mmax;
eq = 0.299^3*M(1)*(i1+1)^4/i1 == 0.255^3*M(1)*i1*(iges/i1+1)^4/(iges/i1); % a(1) == a(2)
sol = solve(eq,i1);
i= [double(sol)];
i = [i iges/i]; % Gear ratios stage 1 and 2

% Empirical values for pinion teeth number from Niemann, Getriebe allgemein, 2003, p.269, case-hardened
i_inter = [1 2 4 8];
z_inter1 = [32 29 25 22]; % high speed
z_inter2 = [26.5 24 20.5 18]; % medium speed

% Interpolation of teeth number for stages 1 and 2 based on empirical values
z(1) = round(interp1(i_inter,z_inter1,i(1),'linear'),0); % high speeds
z(2) = round(interp1(i_inter,z_inter2,i(2),'linear'),0); % medium speeds

% Calculation of gear teeth number(s)
z2 = round(z.*i);


%% Shaft diameter
% Shaft diameter required for seal size and bearing selection

% Equivalent moment due to bending and torsion loading: Factor 1.2 ... 2.5
% Source: Niemann Winter, Maschinenelemente 1, 2019, Formula 17.10
% Mv = 1.5*Mmax; -> used below

% Steel assumption: 16MnCr5 sigma_bW = 520 N/mm^2 from Naunheimer, Fahrzeuggetriebe, 2019, p.478
% Safety factor 1.5, size factor 0.85, surface factor 0.8

sigma_b_zul = 520*0.85*0.8/1.5; % allowable bending stress according to Naunheimer, Fahrzeuggetriebe, 2019, p.478

for j=1:AnzahlStufen+1
    Mv = 1.5*Mmax*prod(i(1:j-1)); % Equivalent moment
    if j > AnzahlStufen % if output shaft
        Mv = Mv/2; % Splitting the moment between left and right
    end
    d(j) = 2.17 * nthroot(Mv*1e3/sigma_b_zul,3); % Shaft diameter according to Naunheimer, Fahrzeuggetriebe, 2019, p.476
end


% Diameter of intermediate shaft according to Niemann, Maschinenelemente 1, 2019, p.494
% Allowable shear stress 190 N/mm^2 based on Taycan shaft determined
dZwischenwelle = 1.72 * nthroot(0.5*Mmax*1e3*iges/190,3);


% Outer diameter of hollow shaft according to Wittel, Roloff/Matek Maschinenelemente, 2021, p.388
syms d1a; % symbolic variable for hollow shaft outer diameter
assume(d1a,'real');
assumeAlso(d1a>0);
eq = 32 * 1.5*Mmax*1e3 / (pi * sigma_b_zul) == (d1a^4 - (dZwischenwelle+1)^4)/d1a;
sol = solve(eq,d1a);
d(1) = double(sol); % Outer diameter of hollow shaft


%% Center distance a, tooth width b, and normal module mn

M(1) = Mmax; % Torque stage 1
M(2) = M(1)*i(1); % Torque stage 2

for j = 1:AnzahlStufen
    % Center distance according to Naunheimer, Fahrzeuggetriebe, 2019, p.288 
    if j==1 % Stage 1
        a(j) = 0.299 * nthroot((M(j)*1e3*(i(j)+1)^4)/i(j),3);
    else % Stage 2
        a(j) = 0.255 * nthroot((M(j)*1e3*(i(j)+1)^4)/i(j),3);
    end
    
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


%% Warning for geometric issues with hollow shaft and pinion on input shaft

if d1(1) - 2*mn(1) < d(1) % Check if base circle diameter < outer diameter of hollow shaft
    disp("- -- -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    disp("Warning: Geometric problem: Hollow shaft too big for calculated pinion in stage 1")
    disp(['Outer diameter hollow shaft: ',sprintf('%0.2f',d(1)),' mm'])
    disp(['Base circle diameter pinion: ',sprintf('%0.2f',d1(1)-2*mn(1)),' mm'])
    disp("- -- -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
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
Gearbox.d1 = d1; % Pitch circle diameter pinion
Gearbox.d2 = d2; % Pitch circle diameter gear(s)
Gearbox.a = a; % Center distance(s)
Gearbox.b = b; % Tooth width(s)
Gearbox.mn = mn; % Normal module(s)
Gearbox.d = d; % Shaft diameter
Gearbox.dLager = dLager; % Bearing inner diameter
Gearbox.DLager = DLager; % Bearing outer diameter
Gearbox.C0Lager = C0Lager; % Static load ratings of bearings
end