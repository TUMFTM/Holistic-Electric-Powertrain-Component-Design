function [P_V] = Getrieberechnung(Gearbox,Mrad,nrad,Temp)
% Function for loss calculation of a 2-stage spur gear with 
% parallel-axis input/output, planar arrangement

%%#codegen

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Oil and lubrication data
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%%% Oil data SAE 70W-75W (API GL4) VW ID3 gearbox oil according to Liqui Moly datasheet
% ny40 = 30.8 mm^2/s and     ny100 = 5.9 mm^2/s    rho15 = 845 kg/m^3
% ny80 = 9.12 mm^2/s  ny90 = 7.25 mm^2/s^2

rho = 820; % kg/m^3: average density of the oil in the operating range

%% Kinematic viscosity of the oil during operation at operating temperature T
% Calculation according to DIN 51563:2011-04
ny40 = 30.8;
T1 = 40;
ny100 = 5.9;
T2 = 100;

W1 = log10(log10(ny40+0.8));
W2 = log10(log10(ny100+0.8));

m = (W1-W2)/(log10(T2+273.15)-log10(T1+273.15));

ny_T = @(T) (10^(10^(m*(log10(T2+273.15)-log10(T+273.15))+W2)))-0.8;

ny = ny_T(Temp);

%% Fixed parameters
XL = 0.75; % Lubrication factor
Qe = 2.5/1000/60; % m^3/s; injection volume flow
v = 0.4; % m/s; injection velocity
Ra = 0.5; % ym: roughness

%%%%%%%%%%%%%%%%%
%% Bearing data
DLager = Gearbox.DLager;
dLager = Gearbox.dLager;
C0Lager = Gearbox.C0Lager;

firstrun = true; % Identifies the first calculation pass, no loss values known yet

%% Preallocation
M = zeros(3,1);
Hv = zeros(2,1);
mym = zeros(2,1);
P_VZP = zeros(2,1);
P_VZ0 = zeros(2,1);
P_VQ = zeros(2,1);
P_VI = zeros(2,1);
P_VZ = zeros(2,1);
P_VD = zeros(2,1);
Mrr = zeros(6,1);
Msl = zeros(6,1);
P_VL = zeros(6,1);

Welle = struct('n',{0 0 0},'P_V',{0 0 0},'L1Fax',{0 0 0},'L2Fax',{0 0 0},'L1Fr',{0 0 0},'L2Fr',{0 0 0},'L1Fy',{0 0 0},'L2Fy',{0 0 0},'L1Fz',{0 0 0},'L2Fz',{0 0 0});
Stufe = struct('P',{0 0},'P_V',{0 0},'Fax',{0 0},'Fu',{0 0},'Frad',{0 0},'dw1',{0 0},'dw2',{0 0});

%% Triple, iterative calculation run (including consideration of
% the loss values calculated in the previous pass)
for k = 1:3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stage 2 (gear side)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha = Gearbox.alpha; % Degrees: Standard pressure angle
    beta = Gearbox.beta(2); % Degrees: Helix angle pitch circle
    mn = Gearbox.mn(2); % mm: Normal module
    z1 = Gearbox.z1(2); % Number of teeth on gear 1
    z2 = Gearbox.z2(2); % Number of teeth on gear 2
    b = Gearbox.b(2); % mm: Tooth width
    

    % Calculate gear data from inputs
    % Script for gear data according to Appendix B of the semester project
    % Sources: Niemann, Getriebe allgemein, 2003, p.58 and 276f.
    %          FZG, Formula collection for the Machine Elements module, 2020, p.102ff.
    %          Gruss, Calculation basics for spur gear teeth (available online) 
        alphawt = atand(tand(alpha)/cosd(beta)); % Degrees: Operating pressure angle
        betab = atand(tand(beta)*cosd(alphawt)); % Degrees: Base helix angle
        mt = mn/cosd(beta); % mm: Transverse module
        
        i_Stufe = z2/z1; % --: Gear ratio; transmission ratio
        
        d1 = mt*z1; % mm: Pitch circle diameter of gear 1
        d2 = mt*z2; % mm: Pitch circle diameter of gear 2
    
        dw1 = d1; % mm: Rolling circle diameter of gear 1
        dw2 = d2; % mm: Rolling circle diameter of gear 2
    
        db1 = d1 * cosd(alphawt); % mm: Base circle diameter of gear 1
        db2 = d2 * cosd(alphawt); % mm: Base circle diameter of gear 2
        
        a = (d1 +d2)/2; % mm: Center distance
        
        haP = mn; % mm: Addendum
        da1 = d1 + 2*haP; % mm: Addendum circle diameter of gear 1
        da2 = d2 + 2*haP; % mm: Addendum circle diameter of gear 2
        
        galpha = 0.5* (sqrt(da1^2 - db1^2) + sqrt(da2^2 - db2^2)) - a*sind(alphawt); % mm: Path of contact
        galphaa1 = 0.5*db1 * (((da1/db1)^2 -1)^0.5 - tand(alphawt)); % mm: Path of contact for pinion addendum
        eps_alpha = galpha / (mt*pi*cosd(alphawt)); % Profile contact ratio
        eps_1 = galphaa1 / (mt*pi*cosd(alphawt)); % Partial contact ratio for pinion addendum
        eps_2 = eps_alpha - eps_1; % Partial contact ratio for gear addendum
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZP (load-dependent gear losses)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     % Calculate speeds and torques in stage 2    
    Welle(3).n = nrad; % Speed of shaft 3 (radial shaft)
    Welle(2).n = Welle(3).n * i_Stufe; % Speed of shaft 2 (after gear stage 2)
    if firstrun % Losses unknown
        M(3) = Mrad; % Torque on shaft 3
    elseif Mrad >= 0 % Driving case (Mrad > 0) -> Torque increases
        M(3) = Mrad + (Stufe(2).P_V + Welle(3).P_V) / (2 * pi * Welle(3).n / 60);
    else % Regeneration case (Mrad < 0) -> Torque decreases
        M(3) = Mrad + Welle(3).P_V / (2 * pi * Welle(3).n / 60);
    end
    
    Stufe(2).P = (Welle(3).n * 2 * pi / 60) * M(3); % Transmitted power in stage 2
    
    % Gear forces
    Fu = 2 * M(3) / (d2 * 1e-3); % Circumferential force
    Fax = Fu * tand(beta); % Axial force
    Frad = abs(Fu) * tand(alphawt); % Radial force
    
    % Intermediate values for P_VZP calculation
    Fbt = abs(Fu / cosd(alphawt)); % Circumferential force normal to the tooth
    vt = Welle(3).n / 60 * pi * d2 * 1e-3; % Circumferential velocity of the gear
    vsigC = 2 * vt * sind(alphawt); % Summation velocity at rolling contact point
    etaOil = ny * rho / 1e3; % Dynamic viscosity of the oil
    rhoredC = 0.5 * d1 * sind(alphawt) * i_Stufe(1) / (cosd(betab) * (i_Stufe(1) + 1)); % Equivalent curvature radius at rolling contact point
    
    % Gear loss factor Hv according to Ohlendorf
    Hv(2) = pi * (i_Stufe(1) + 1) / (z1 * i_Stufe(1) * cosd(betab)) * (1 - eps_alpha + eps_1^2 + eps_2^2);
    
    % Mean gear friction coefficient according to Schlenk
    mym(2) = 0.048 * (Fbt / b / (vsigC * rhoredC))^0.2 * etaOil^-0.05 * XL * Ra^0.25;
    
    % Load-dependent power losses in stage 2
    P_VZP(2) = abs(Stufe(2).P) * mym(2) * Hv(2);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZ0 (load-independent gear losses) with injection lubrication
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Squeeze losses according to Mauz
    T_Q = 0; % No squeeze losses if injected at the outlet (assumption from WTplus)
    P_VQ(2) = T_Q * Welle(2).n * 2 * pi / 60; % = 0 because T_Q = 0
    
    % Impulse losses according to Ariura
    C1 = 1; 
    P_VI(2) = d2 / 2 * 1e-3 * rho * Qe * (vt + v) * C1 * Welle(3).n * 2 * pi / 60;
    
    % Ventilation losses
    % P_VV = 1.37e-9 * vt^1.9 * d1^1.6 * b^0.52 * mn^0.69 * na * 2 * pi / 60
    % Neglected
    
    % Load-independent gear losses
    P_VZ0(2) = P_VQ(2) + P_VI(2);

    % Total gear losses
    P_VZ(2) = P_VZP(2) + P_VZ0(2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Bearing forces on shaft 3 (see force diagram in appendix of SA)
    a = 80 + b / 2; % Distance adjustment for differential
    b1 = 0; 
    c = 20 + b / 2; 
    Welle(3).L1Fax = Fax;
    Welle(3).L2Fax = 0;
    Welle(3).L1Fy = -(Frad * c + Fax * dw2 / 2) / (a + b1 + c);
    Welle(3).L2Fy = -(Frad + Welle(3).L1Fy);
    Welle(3).L2Fz = (Fu * a) / (a + b1 + c);
    Welle(3).L1Fz = (Fu - Welle(3).L2Fz);
    Welle(3).L1Fr = sqrt(Welle(3).L1Fy^2 + Welle(3).L1Fz^2);
    Welle(3).L2Fr = sqrt(Welle(3).L2Fy^2 + Welle(3).L2Fz^2);
    
    % Save stage data
    Stufe(2).Fax = Fax;
    Stufe(2).Fu = Fu;
    Stufe(2).Frad = Frad;
    Stufe(2).dw1 = dw1;
    
    % Calculate torque on shaft 2
    if firstrun % Losses unknown
        M(2) = M(3) / i_Stufe;
    elseif Mrad >= 0 % Driving case, losses increase power
        M(2) = M(3) / i_Stufe + (Stufe(1).P_V + Welle(2).P_V) / (2 * pi * Welle(2).n / 60);
    else % Regeneration case, losses reduce negative power
        M(2) = M(3) / i_Stufe + (Stufe(2).P_V + Welle(2).P_V) / (2 * pi * Welle(2).n / 60);
    end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stage 1 (motor side)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = Gearbox.beta(1); % Degrees: Helix angle pitch circle
    mn = Gearbox.mn(1); % mm: Normal module
    z1 = Gearbox.z1(1); % Number of teeth on gear 1
    z2 = Gearbox.z2(1); % Number of teeth on gear 2
    b = Gearbox.b(1); % mm: Tooth width
    
    % Calculate gear data from inputs
    % Script for gear data according to Appendix B of the semester project
    % Sources: Niemann, Getriebe allgemein, 2003, p.58 and 276f.
    %          FZG, Formula collection for the Machine Elements module, 2020, p.102ff.
    %          Gruss, Calculation basics for spur gear teeth (available online) 
        alphawt = atand(tand(alpha)/cosd(beta)); % Degrees: Operating pressure angle
        betab = atand(tand(beta)*cosd(alphawt)); % Degrees: Base helix angle
        mt = mn/cosd(beta); % mm: Transverse module
        
        i_Stufe = z2/z1; % --: Gear ratio; transmission ratio
        
        d1 = mt*z1; % mm: Pitch circle diameter of gear 1
        d2 = mt*z2; % mm: Pitch circle diameter of gear 2
    
        dw1 = d1; % mm: Rolling circle diameter of gear 1
        dw2 = d2; % mm: Rolling circle diameter of gear 2
    
        db1 = d1 * cosd(alphawt); % mm: Base circle diameter of gear 1
        db2 = d2 * cosd(alphawt); % mm: Base circle diameter of gear 2
        
        a = (d1 + d2)/2; % mm: Center distance
        
        haP = mn; % mm: Addendum
        da1 = d1 + 2*haP; % mm: Addendum circle diameter of gear 1
        da2 = d2 + 2*haP; % mm: Addendum circle diameter of gear 2
        
        galpha = 0.5 * (sqrt(da1^2 - db1^2) + sqrt(da2^2 - db2^2)) - a * sind(alphawt); % mm: Path of contact
        galphaa1 = 0.5 * db1 * (((da1/db1)^2 - 1)^0.5 - tand(alphawt)); % mm: Path of contact for pinion addendum
        eps_alpha = galpha / (mt * pi * cosd(alphawt)); % Profile contact ratio
        eps_1 = galphaa1 / (mt * pi * cosd(alphawt)); % Partial contact ratio for pinion addendum
        eps_2 = eps_alpha - eps_1; % Partial contact ratio for gear addendum
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZP (load-dependent gear losses) in stage 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Speed of shaft 1 and transmitted power in stage 1
    Welle(1).n = Welle(2).n * i_Stufe;
    Stufe(1).P = (Welle(2).n * 2 * pi / 60) * M(2);
    
    % Gear forces
    Fu = 2 * M(2) / (d2 * 1e-3); % Circumferential force
    Fax = -Fu * tand(beta); % Axial force
    Frad = abs(Fu) * tand(alphawt); % Radial force
    
    % Intermediate values for P_VZP calculation
    Fbt = abs(Fu / cosd(alphawt)); % Circumferential force normal to the tooth
    vt = Welle(2).n / 60 * pi * d2 * 1e-3; % Circumferential velocity of the gear
    vsigC = 2 * vt * sind(alphawt); % Summation velocity at rolling contact point
    etaOil = rho * ny / 1e3; % Dynamic viscosity of the oil
    rhoredC = 0.5 * d1 * sind(alphawt) * i_Stufe(1) / (cosd(betab) * (i_Stufe(1) + 1)); % Equivalent curvature radius at rolling contact point
    
    % Gear loss factor Hv according to Ohlendorf
    Hv(1) = pi * (i_Stufe(1) + 1) / (z1 * i_Stufe(1) * cosd(betab)) * (1 - eps_alpha + eps_1^2 + eps_2^2);
    
    % Mean gear friction coefficient according to Schlenk
    mym(1) = 0.048 * (Fbt / b / (vsigC * rhoredC))^0.2 * etaOil^-0.05 * XL * Ra^0.25;
    
    % Load-dependent power losses in stage 1
    P_VZP(1) = abs(Stufe(1).P) * mym(1) * Hv(1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZ0 (load-independent gear losses) with injection lubrication
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % Squeeze losses according to Mauz
    T_Q = 0; % No squeeze losses if injected at the outlet (WTplus assumption)
    P_VQ(1) = T_Q * Welle(1).n * 2 * pi / 60;
    
    % Impulse losses according to Ariura
    C1 = 0.85; 
    P_VI(1) = d2 / 2 * 1e-3 * rho * Qe * (vt + v) * C1 * Welle(2).n * 2 * pi / 60;

    % Ventilation losses P_VV neglected
    
    % Load-independent gear losses in stage 1:
    P_VZ0(1) = P_VQ(1) + P_VI(1);

    % Total gear losses in stage 1:
    P_VZ(1) = P_VZP(1) + P_VZ0(1);
    
    % Save stage data
    Stufe(1).Fax = Fax;
    Stufe(1).Fu = Fu;
    Stufe(1).Frad = Frad;
    Stufe(1).dw2 = dw2;
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Bearing forces on shaft 2 (see force diagram in the appendix of SA)
    a = 8 + Gearbox.b(1)/2; 
    b1 = Gearbox.b(1)/2 + 5 + Gearbox.b(2)/2; 
    c = Gearbox.b(2)/2 + 8;
    Welle(2).L1Fax = - Stufe(1).Fax - Stufe(2).Fax;
    Welle(2).L2Fax = 0;
    Welle(2).L1Fy = (Stufe(1).Fax * Stufe(1).dw2 / 2 - Stufe(1).Frad * (b1 + c) + Stufe(2).Frad * c - Stufe(2).Fax * Stufe(2).dw1 / 2) / (a + b1 + c);
    Welle(2).L2Fy = -(Stufe(1).Frad - Stufe(2).Frad + Welle(2).L1Fy);
    Welle(2).L2Fz = -(Stufe(1).Fu * a + Stufe(2).Fu * (a + b1)) / (a + b1 + c);
    Welle(2).L1Fz = -(Stufe(1).Fu + Stufe(2).Fu + Welle(2).L2Fz);
    Welle(2).L1Fr = sqrt(Welle(2).L1Fy^2 + Welle(2).L1Fz^2);
    Welle(2).L2Fr = sqrt(Welle(2).L2Fy^2 + Welle(2).L2Fz^2);
    
    
    %% Bearing forces on shaft 1 (see force diagram in the appendix of SA)
    a = 8 + b / 2; 
    b1 = 0; 
    c = b / 2 + 8; 
    Welle(1).L1Fax = Fax;
    Welle(1).L2Fax = 0;
    Welle(1).L1Fy = (Fax * dw1 / 2 + Frad * c) / (a + b1 + c);
    Welle(1).L2Fy = Frad - Welle(1).L1Fy;
    Welle(1).L2Fz = Fu * a / (a + b1 + c);
    Welle(1).L1Fz = Fu - Welle(1).L2Fz;
    Welle(1).L1Fr = sqrt(Welle(1).L1Fy^2 + Welle(1).L1Fz^2);
    Welle(1).L2Fr = sqrt(Welle(1).L2Fy^2 + Welle(1).L2Fz^2);
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VD Seal losses according to ISO 14179-2 (Germany)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Seal of motor input shaft
    P_VD(1) = 7.69e-6 * (Gearbox.d(1)^2 * Welle(1).n);
    
    % Seals (2x) for output shaft (wheel side)
    P_VD(2) = 2 * 7.69e-6 * (Gearbox.d(3)^2 * Welle(3).n);

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VL Bearing losses SKF2020 (hardcoded for 160/161 deep groove ball bearings)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Formulas according to SKF, The SKF model for calculating the frictional moment
    % Online: https://cdn.skfmediahub.skf.com/api/public/0901d196809bc183/pdf_preview_medium/0901d196809bc183_pdf_preview_medium.pdf

    % Bearing speeds = Shaft speeds
    nLager = [Welle(1).n Welle(2).n Welle(3).n];
    
    for i = 1:3 % Iterate over shafts
        for j = 1:2 % Iterate over bearings per shaft
            index = 2 * (i - 1) + j;

            % Read specific bearing values from design
            DLag = DLager(i); % Outer diameter of the bearing
            dLag = dLager(i); % Inner diameter of the bearing
            C0 = C0Lager(i); % Static load rating of the bearing
            nLag = nLager(i); % Bearing speed
            
            dm = (DLag + dLag) / 2; % Average bearing diameter
            
            % Assign acting forces
            if j == 1 % Bearing 1
                Fr = Welle(i).L1Fr; 
                Fa = abs(Welle(i).L1Fax);
            else % Bearing 2
                Fr = Welle(i).L2Fr; 
                Fa = abs(Welle(i).L2Fax);
            end
            
            % Coefficients for 160/161 series deep groove ball bearings
            S1 = 4.63e-3; 
            S2 = 4.25;
            R1 = 4.3e-7; 
            R2 = 1.7;
            
            % Lubrication film thickness factor psi_ish
            psiish = 1 / (1 + 1.84e-9 * (nLag * dm)^1.28 * ny^0.64);

            % Lubrication displacement factor psi_rs
            psirs = 1 / exp(6e-8 * ny * nLag * (DLag + dLag) * sqrt(3.1 / (2 * (DLag - dLag))));

            % Base rolling friction G_rr, Base sliding friction G_sl
            if Fa == 0 % Without axial load
                Grr = R1 * dm^1.96 * Fr^0.54;
                Gsl = S1 * dm^(-0.26) * Fr^(5/3);
            else % With axial load
                alphaf = 24.6 * (Fa / C0)^0.24;
                Grr = R1 * dm^1.96 * (Fr + R2 / sind(alphaf) * Fa)^0.54;
                Gsl = S1 * dm^-0.145 * (Fr^5 + S2 * dm^1.5 / sind(alphaf) * Fa^4)^(1/3);
            end

            % Rolling friction moment M_rr
            Mrr(index) = psiish * psirs * Grr * (ny * nLag)^0.6;
            
            % Weighting factor for sliding friction coefficient psi_bl
            psibl = 1 / exp(2.6e-8 * (nLag * ny)^1.4 * dm);

            % Sliding friction coefficient my_sl
            mysl = psibl * 0.12 + (1 - psibl) * 0.04; % SKF online calculator uses 0.05 instead of 0.04

            % Sliding friction moment M_sl
            Msl(index) = Gsl * mysl;
            
            % Bearing power losses
            P_VL(index) = (Mrr(index) + Msl(index)) * 1e-3 * nLag * 2 * pi / 60;
        end
    end
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%
    % Assignment of losses to the corresponding shaft or stage
    % Necessary to consider them at the correct point in the next iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Stufe(1).P_V = P_VZ(1); % Assign gear losses to stage 1
    Stufe(2).P_V = P_VZ(2); % Assign gear losses to stage 2
    Welle(1).P_V = P_VD(1) + sum(P_VL(1:2)); % Assign total losses to shaft 1
    Welle(2).P_V = sum(P_VL(3:4)); % Assign bearing losses to shaft 2
    Welle(3).P_V = P_VD(2) + sum(P_VL(5:6)); % Assign total losses to shaft 3
    
    firstrun = false; % In subsequent iterations, losses are known
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Total losses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_V = sum(P_VD) + sum(P_VZ) + sum(P_VL); % Calculate total loss

% To display the loss components in the Command Window, activate the following:
% disp(['Total gear loss power: ',sprintf('%0.2f',sum(P_VZ)),' W'])
% disp(['Total seal loss power: ',sprintf('%0.2f',sum(P_VD)),' W'])
% disp(['Total bearing loss power: ',sprintf('%0.2f',sum(P_VL)),' W'])
% disp(['Total gearbox loss power: ',sprintf('%0.2f',P_V),' W'])

end