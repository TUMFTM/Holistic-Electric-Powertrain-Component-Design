function [P_V] = Getrieberechnung_koaxial(Gearbox,Mrad,nrad,Temp)
% Function for calculating the losses of a 2-stage spur gear with 
% coaxial input/output and planar arrangement.

%% Oil and lubrication data

%%%%%%%%%%%%%%%%
%%% Oil data SAE 70W-75W (API GL4) VW ID3 gear oil
% ny40 = 30.8 mm^2/s and ny100 = 5.9 mm^2/s, rho15 = 845 kg/m^3
% ny80 = 9.12 mm^2/s, ny90 = 7.25 mm^2/s

rho = 820; % kg/m^3: Average density of the oil in the operating range.

%% Kinematic viscosity of the oil in operation at operating temperature T
% Calculation according to DIN 51563:2011-04
ny40 = 30.8; % Kinematic viscosity at 40°C
T1 = 40;     % Reference temperature 1 in °C
ny100 = 5.9; % Kinematic viscosity at 100°C
T2 = 100;    % Reference temperature 2 in °C

W1 = log10(log10(ny40 + 0.8));
W2 = log10(log10(ny100 + 0.8));

m = (W1 - W2) / (log10(T2 + 273.15) - log10(T1 + 273.15));

% Viscosity as a function of temperature
ny_T = @(T) (10^(10^(m * (log10(T2 + 273.15) - log10(T + 273.15)) + W2))) - 0.8;

ny = ny_T(Temp);

%% Fixed parameters
XL = 0.75; % Lubrication factor
Qe = 2.5 / 1000 / 60; % m^3/s: Injection flow rate
v = 0.4; % m/s: Injection speed
Ra = 0.5; % µm: Roughness

%% Bearing data
DLager = Gearbox.DLager; % Outer diameter of bearings
dLager = Gearbox.dLager; % Inner diameter of bearings
C0Lager = Gearbox.C0Lager; % Static load rating of bearings

firstrun = true; % Identifies the first calculation run, losses are not yet known.

%% Preallocation of variables
M = zeros(3, 1);
Hv = zeros(2, 1);
mym = zeros(2, 1);
P_VZP = zeros(2, 1);
P_VZ0 = zeros(2, 1);
P_VQ = zeros(2, 1);
P_VI = zeros(2, 1);
P_VZ = zeros(2, 1);
P_VD = zeros(2, 1);
Mrr = zeros(6, 1);
Msl = zeros(6, 1);
P_VL = zeros(6, 1);

% Preallocating structures for shafts and stages
Welle = struct('n', {0 0 0}, 'P_V', {0 0 0}, 'L1Fax', {0 0 0}, 'L2Fax', {0 0 0}, ...
    'L1Fr', {0 0 0}, 'L2Fr', {0 0 0}, 'L1Fy', {0 0 0}, 'L2Fy', {0 0 0}, ...
    'L1Fz', {0 0 0}, 'L2Fz', {0 0 0});

Stufe = struct('P', {0 0}, 'P_V', {0 0}, 'Fax', {0 0}, 'Fu', {0 0}, ...
    'Frad', {0 0}, 'dw1', {0 0}, 'dw2', {0 0});

%% Triple iterative calculation loop (including consideration of 
% loss values calculated in the previous iteration).
for k = 1:3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stage 2 (wheel side)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha = Gearbox.alpha; % Degrees: Normal pressure angle
    beta = Gearbox.beta(2); % Degrees: Helix angle (pitch circle)
    mn = Gearbox.mn(2); % mm: Normal module
    z1 = Gearbox.z1(2); % Number of teeth on gear 1
    z2 = Gearbox.z2(2); % Number of teeth on gear 2
    b = Gearbox.b(2); % mm: Face width

    % Calculate gear data from inputs
    % Script for gear data based on Appendix B of the semester paper
    % Sources: Niemann, Gears in general, 2003, pp. 58 and 276f.
    %          FZG, Formula collection for the machine elements module, 2020, pp. 102ff.
    %          Gruss, Fundamentals of spur gear calculations (available online)
    alphawt = atand(tand(alpha) / cosd(beta)); % Degrees: Operating pressure angle
    betab = atand(tand(beta) * cosd(alphawt)); % Degrees: Base helix angle
    mt = mn / cosd(beta); % mm: Transverse module
    
    i_Stufe = z2 / z1; % Gear ratio of the stage

    d1 = mt * z1; % mm: Pitch diameter of gear 1
    d2 = mt * z2; % mm: Pitch diameter of gear 2

    dw1 = d1; % mm: Rolling circle diameter gear 1
    dw2 = d2; % mm: Rolling circle diameter gear 2

    db1 = d1 * cosd(alphawt); % mm: Base diameter gear 1
    db2 = d2 * cosd(alphawt); % mm: Base diameter gear 2

    a = (d1 + d2) / 2; % mm: Center distance

    haP = mn; % mm: Addendum
    da1 = d1 + 2 * haP; % mm: Tip diameter gear 1
    da2 = d2 + 2 * haP; % mm: Tip diameter gear 2

    galpha = 0.5 * (sqrt(da1^2 - db1^2) + sqrt(da2^2 - db2^2)) - a * sind(alphawt); % mm: Contact line length
    galphaa1 = 0.5 * db1 * (((da1 / db1)^2 - 1)^0.5 - tand(alphawt)); % mm: Contact line length (pinion addendum)
    eps_alpha = galpha / (mt * pi * cosd(alphawt)); % Profile contact ratio
    eps_1 = galphaa1 / (mt * pi * cosd(alphawt)); % Partial contact ratio (pinion)
    eps_2 = eps_alpha - eps_1; % Partial contact ratio (wheel)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZP (Load-dependent gear losses)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate rotational speeds and torques in stage 2
    Welle(3).n = nrad; % Speed of shaft 3
    Welle(2).n = Welle(3).n * i_Stufe; % Speed of shaft 2 based on gear ratio
    if firstrun % Losses are unknown
        M(3) = Mrad; % Torque on shaft 3
    elseif Mrad >= 0 % Driving mode (Mrad > 0) -> Torque increases
        M(3) = Mrad + (Stufe(2).P_V + Welle(3).P_V) / (2 * pi * Welle(3).n / 60);
    else % Regeneration mode (Mrad < 0) -> Torque decreases
        M(3) = Mrad + Welle(3).P_V / (2 * pi * Welle(3).n / 60);
    end
    
    Stufe(2).P = (Welle(3).n * 2 * pi / 60) * M(3); % Power transmitted in stage 2
    
    % Gear forces
    Fu = 2 * M(3) / (d2 * 1e-3); % Tangential force
    Fax = Fu * tand(beta); % Axial force
    Frad = abs(Fu) * tand(alphawt); % Radial force
    
    % Intermediate values for P_VZP calculation
    Fbt = abs(Fu / cosd(alphawt)); % Tangential force normal to the tooth
    vt = Welle(3).n / 60 * pi * d2 * 1e-3; % Tangential velocity of the gear
    vsigC = 2 * vt * sind(alphawt); % Combined velocity at the rolling point
    etaOil = rho * ny / 1e3; % Dynamic viscosity of the oil
    rhoredC = 0.5 * d1 * sind(alphawt) * i_Stufe(1) / (cosd(betab) * (i_Stufe(1) + 1)); % Equivalent curvature radius at the rolling point
    
    % Gear loss factor Hv according to Ohlendorf
    Hv(2) = pi * (i_Stufe(1) + 1) / (z1 * i_Stufe(1) * cosd(betab)) * (1 - eps_alpha + eps_1^2 + eps_2^2);
    
    % Average friction coefficient according to Schlenk
    mym(2) = 0.048 * (Fbt / b / (vsigC * rhoredC))^0.2 * etaOil^-0.05 * XL * Ra^0.25;
    
    % Load-dependent power loss in stage 2
    P_VZP(2) = abs(Stufe(2).P) * mym(2) * Hv(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZ0 (Load-independent gear losses) with injection lubrication
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Squeeze losses according to Mauz
    T_Q = 0; % No squeeze losses if lubricated at the outlet (simplified assumption)
    P_VQ(2) = T_Q * Welle(2).n * 2 * pi / 60;
    
    % Impulse losses according to Ariura
    C1 = 1;
    P_VI(2) = d2 / 2 * 1e-3 * rho * Qe * (vt + v) * C1 * Welle(3).n * 2 * pi / 60;
    
    % Ventilation losses (neglected)
    % P_VV = 1.37e-9 * vt^1.9 * d1^1.6 * b^0.52 * mn^0.69 * na * 2 * pi / 60;
    
    % Load-independent gear losses
    P_VZ0(2) = P_VQ(2) + P_VI(2);
    
    % Total gear losses
    P_VZ(2) = P_VZP(2) + P_VZ0(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Bearing forces on shaft 3 (See force diagram in the appendix)
    a = 70 + b / 2; % Distance a
    b1 = 0; % Distance b1
    c = 13 + b / 2; % Distance c for differential clearance
    Welle(3).L1Fax = Fax;
    Welle(3).L2Fax = 0;
    Welle(3).L1Fy = (Frad * c + Fax * dw2 / 2) / (a + b1 + c);
    Welle(3).L2Fy = (Frad - Welle(3).L1Fy);
    Welle(3).L2Fz = -(Fu * a) / (a + b1 + c);
    Welle(3).L1Fz = (-Fu - Welle(3).L2Fz);
    Welle(3).L1Fr = sqrt(Welle(3).L1Fy^2 + Welle(3).L1Fz^2);
    Welle(3).L2Fr = sqrt(Welle(3).L2Fy^2 + Welle(3).L2Fz^2);
    
    % Store stage 2 data
    Stufe(2).Fax = Fax;
    Stufe(2).Fu = Fu;
    Stufe(2).Frad = Frad;
    Stufe(2).dw1 = dw1;
    
    % Torque calculation on shaft 2
    if firstrun % Losses unknown
        M(2) = M(3) / i_Stufe;
    elseif Mrad >= 0 % Driving mode, losses increase power
        M(2) = M(3) / i_Stufe + (Stufe(1).P_V + Welle(2).P_V) / (2 * pi * Welle(2).n / 60);
    else % Regeneration mode, losses decrease negative power
        M(2) = M(3) / i_Stufe + (Stufe(2).P_V + Welle(2).P_V) / (2 * pi * Welle(2).n / 60);
    end
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stage 1 (motor side)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = Gearbox.beta(1); % Degrees: Helix angle at the pitch circle
    mn = Gearbox.mn(1); % mm: Normal module
    z1 = Gearbox.z1(1); % Number of teeth on gear 1
    z2 = Gearbox.z2(1); % Number of teeth on gear 2
    b = Gearbox.b(1); % mm: Tooth width
    
    % Calculate gear data from inputs
    % Script for gear data according to Appendix B of the semester project
    % Sources: Niemann, General Gear Design, 2003, pp. 58 and 276ff.
    %          FZG, Formula Collection for the Machine Elements Module, 2020, pp. 102ff.
    %          Gruss, Calculation Basics of Spur Gear Teeth (available online) 
    alphawt = atand(tand(alpha) / cosd(beta)); % Degrees: Operating pressure angle
    betab = atand(tand(beta) * cosd(alphawt)); % Degrees: Helix angle at the base circle
    mt = mn / cosd(beta); % mm: Transverse module
    
    i_Stufe = z2 / z1; % --: Gear ratio; reduction ratio
    
    d1 = mt * z1; % mm: Pitch diameter of gear 1
    d2 = mt * z2; % mm: Pitch diameter of gear 2
    
    dw1 = d1; % mm: Rolling circle diameter of gear 1
    dw2 = d2; % mm: Rolling circle diameter of gear 2
    
    db1 = d1 * cosd(alphawt); % mm: Base circle diameter of gear 1
    db2 = d2 * cosd(alphawt); % mm: Base circle diameter of gear 2
    
    a = (d1 + d2) / 2; % mm: Center distance
    
    haP = mn; % mm: Addendum height
    da1 = d1 + 2 * haP; % mm: Addendum circle diameter of gear 1
    da2 = d2 + 2 * haP; % mm: Addendum circle diameter of gear 2
    
    galpha = 0.5 * (sqrt(da1^2 - db1^2) + sqrt(da2^2 - db2^2)) - a * sind(alphawt); % mm: Contact line length
    galphaa1 = 0.5 * db1 * (((da1 / db1)^2 - 1)^0.5 - tand(alphawt)); % mm: Contact line length at pinion addendum
    eps_alpha = galpha / (mt * pi * cosd(alphawt)); % Contact ratio
    eps_1 = galphaa1 / (mt * pi * cosd(alphawt)); % Partial contact ratio for pinion addendum
    eps_2 = eps_alpha - eps_1; % Partial contact ratio for gear addendum
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZP (Load-dependent gear losses) in Stage 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Shaft 1 speed and transmitted power in Stage 1
    Welle(1).n = Welle(2).n * i_Stufe; % Shaft 1 speed
    Stufe(1).P = (Welle(2).n * 2 * pi / 60) * M(2); % Power transmitted in Stage 1
    
    % Gear forces
    Fu = 2 * M(2) / (d2 * 1e-3); % Tangential force
    Fax = Fu * tand(beta); % Axial force
    Frad = abs(Fu) * tand(alphawt); % Radial force
    
    % Intermediate values for P_VZP calculation
    Fbt = abs(Fu / cosd(alphawt)); % Tangential force normal to the tooth
    vt = Welle(2).n / 60 * pi * d2 * 1e-3; % Tangential velocity of the gear
    vsigC = 2 * vt * sind(alphawt); % Combined velocity at the rolling point
    etaOil = rho * ny / 1e3; % Dynamic viscosity of the oil
    rhoredC = 0.5 * d1 * sind(alphawt) * i_Stufe(1) / (cosd(betab) * (i_Stufe(1) + 1)); % Equivalent curvature radius at the rolling point
    
    % Gear loss factor Hv according to Ohlendorf
    Hv(1) = pi * (i_Stufe(1) + 1) / (z1 * i_Stufe(1) * cosd(betab)) * (1 - eps_alpha + eps_1^2 + eps_2^2);
    
    % Average friction coefficient according to Schlenk
    mym(1) = 0.048 * (Fbt / b / (vsigC * rhoredC))^0.2 * etaOil^-0.05 * XL * Ra^0.25;
    
    % Load-dependent power loss in Stage 1
    P_VZP(1) = abs(Stufe(1).P) * mym(1) * Hv(1);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZ0 (Load-independent gear losses) with injection lubrication
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compression losses according to Mauz
    T_Q = 0; % No compression losses when injected into the outlet (simplified assumption)
    P_VQ(1) = T_Q * Welle(1).n * 2 * pi / 60;
    
    % Impulse losses according to Ariura
    C1 = 0.85;
    P_VI(1) = d2 / 2 * 1e-3 * rho * Qe * (vt + v) * C1 * Welle(2).n * 2 * pi / 60;
    
    % Ventilation losses P_VV are neglected
    
    % Load-independent gear losses for Stage 1:
    P_VZ0(1) = P_VQ(1) + P_VI(1);
    
    % Total gear losses for Stage 1:
    P_VZ(1) = P_VZP(1) + P_VZ0(1);
    
    % Save stage data
    Stufe(1).Fax = Fax;
    Stufe(1).Fu = Fu;
    Stufe(1).Frad = Frad;
    Stufe(1).dw2 = dw2;
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bearing forces on Shaft 2 (see force diagram in the appendix of the thesis)
    a = 8+Gearbox.b(2)/2; b1 = Gearbox.b(2)/2+29+b/2; c = b/2+8;
    Welle(2).L1Fax = Stufe(1).Fax - Stufe(2).Fax;
    Welle(2).L2Fax = 0;
    Welle(2).L1Fy = (-Stufe(1).Fax*Stufe(1).dw2/2 - Stufe(1).Frad*(c) - Stufe(2).Frad*(b1+c) + Stufe(2).Fax*Stufe(2).dw1/2)/(a+b1+c);
    Welle(2).L2Fy = (-Stufe(1).Frad - Stufe(2).Frad - Welle(2).L1Fy);
    Welle(2).L2Fz = (-Stufe(1).Fu*(a+b1) + Stufe(2).Fu*(a)) / (a+b1+c);
    Welle(2).L1Fz = (-Stufe(1).Fu + Stufe(2).Fu - Welle(2).L2Fz);
    Welle(2).L1Fr = sqrt(Welle(2).L1Fy^2 + Welle(2).L1Fz^2);
    Welle(2).L2Fr = sqrt(Welle(2).L2Fy^2 + Welle(2).L2Fz^2);
    
    
    %% Bearing forces on Shaft 1 (see force diagram in the appendix of the thesis)
    a = 8+b/2; b1 = 0; c = b/2+8; 
    Welle(1).L1Fax = 0;
    Welle(1).L2Fax = -Fax;
    Welle(1).L1Fy = (-Fax*dw1/2 + Frad*c)/(a+b1+c);
    Welle(1).L2Fy = Frad - Welle(1).L1Fy;
    Welle(1).L2Fz = Fu*a/(a+b1+c);
    Welle(1).L1Fz = Fu- Welle(1).L2Fz;
    Welle(1).L1Fr = sqrt(Welle(1).L1Fy^2 + Welle(1).L1Fz^2);
    Welle(1).L2Fr = sqrt(Welle(1).L2Fy^2 + Welle(1).L2Fz^2);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VD Seal Losses according to ISO 14179-2 (Germany)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Seal losses for the motor's input hollow shaft
    P_VD(1) = 7.69e-6 * (Gearbox.d(1)^2 * Welle(1).n);
    
    % Seal losses for the output shaft (2x seals for each wheel axle)
    P_VD(2) = 2 * 7.69e-6 * (Gearbox.d(3)^2 * Welle(3).n);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VL Bearing Losses SKF2020 (Hardcoded for 160/161 Ball Bearings)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Formulas based on SKF, The SKF model for calculating frictional moments
    % Source: https://cdn.skfmediahub.skf.com/api/public/0901d196809bc183/pdf_preview_medium/0901d196809bc183_pdf_preview_medium.pdf
    
    % Bearing speeds are equivalent to shaft speeds
    nLager = [Welle(1).n, Welle(2).n, Welle(3).n];
    
    for i = 1:3 % Iterate over the shafts
        for j = 1:2 % Iterate over the two bearings per shaft
            index = 2 * (i-1) + j; % Calculate the unique index for each bearing
    
            % Retrieve specific bearing parameters from the design
            DLag = DLager(i); % Outer diameter of the bearing
            dLag = dLager(i); % Inner diameter of the bearing
            C0 = C0Lager(i); % Static load rating of the bearing
            nLag = nLager(i); % Rotational speed of the bearing
    
            dm = (DLag + dLag) / 2; % Average bearing diameter
    
            % Assign acting forces
            if j == 1 % For bearing 1
                Fr = Welle(i).L1Fr; % Radial force
                Fa = abs(Welle(i).L1Fax); % Axial force
            else % For bearing 2
                Fr = Welle(i).L2Fr; % Radial force
                Fa = abs(Welle(i).L2Fax); % Axial force
            end
    
            % Coefficients for 160/161 ball bearing series
            S1 = 4.63e-3;
            S2 = 4.25;
            R1 = 4.3e-7;
            R2 = 1.7;
    
            % Lubricant film thickness factor `psi_ish`
            psiish = 1 / (1 + 1.84e-9 * (nLag * dm)^1.28 * ny^0.64);
    
            % Lubricant displacement factor `psi_rs`
            psirs = 1 / exp(6e-8 * ny * nLag * (DLag + dLag) * sqrt(3.1 / (2 * (DLag - dLag))));
    
            % Base rolling friction and sliding friction values
            if Fa == 0 % Without axial load
                Grr = R1 * dm^1.96 * Fr^0.54; % Rolling friction
                Gsl = S1 * dm^(-0.26) * Fr^(5/3); % Sliding friction
            else % With axial load
                alphaf = 24.6 * (Fa / C0)^0.24; % Load angle factor
                Grr = R1 * dm^1.96 * (Fr + R2 / sind(alphaf) * Fa)^0.54; % Rolling friction
                Gsl = S1 * dm^-0.145 * (Fr^5 + S2 * dm^1.5 / sind(alphaf) * Fa^4)^(1/3); % Sliding friction
            end
    
            % Rolling friction moment `M_rr`
            Mrr(index) = psiish * psirs * Grr * (ny * nLag)^0.6;
    
            % Sliding friction coefficient weighting factor `psi_bl`
            psibl = 1 / exp(2.6e-8 * (nLag * ny)^1.4 * dm);
    
            % Sliding friction coefficient `my_sl`
            mysl = psibl * 0.12 + (1 - psibl) * 0.04; % SKF online calculator uses 0.05 instead of 0.04
    
            % Sliding friction moment `M_sl`
            Msl(index) = Gsl * mysl;
    
            % Bearing loss power `P_VL`
            P_VL(index) = (Mrr(index) + Msl(index)) * 1e-3 * nLag * 2 * pi / 60;
        end
    end



   %%%%%%%%%%%%%%%%%%%%%%%%%
    % Allocation of losses to the corresponding shaft or stage
    % Necessary to account for them correctly in the next iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Stufe(1).P_V = P_VZ(1);
    Stufe(2).P_V = P_VZ(2);
    Welle(1).P_V = P_VD(1) + sum(P_VL(1:2));
    Welle(2).P_V = sum(P_VL(3:4));
    Welle(3).P_V = P_VD(2) + sum(P_VL(5:6));
    
    firstrun = false; % In subsequent iterations, the losses are known
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Total Loss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_V = sum(P_VD) + sum(P_VZ) + sum(P_VL);

% To display the loss components in the Command Window, activate:
% disp(['Total Gear Loss Power: ', sprintf('%0.2f', sum(P_VZ)), ' W'])
% disp(['Total Seal Loss Power: ', sprintf('%0.2f', sum(P_VD)), ' W'])
% disp(['Total Bearing Loss Power: ', sprintf('%0.2f', sum(P_VL)), ' W'])
% disp(['Total Gearbox Loss Power: ', sprintf('%0.2f', P_V), ' W'])

end