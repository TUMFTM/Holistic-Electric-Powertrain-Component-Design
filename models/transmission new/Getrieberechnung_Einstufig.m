function [P_V] = Getrieberechnung_Einstufig(Gearbox,Mrad,nrad,Temp)
% Function for loss calculation of a single-stage spur gear with
% coaxial input/output, planar arrangement

%%#codegen

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Oil and lubrication data
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%%% Oil data SAE 70W-75W (API GL4) VW ID3 gearbox oil according to Liqui Moly datasheet
% ny40 = 30.8 mm^2/s and ny100 = 5.9 mm^2/s, rho15 = 845 kg/m^3
% ny80 = 9.12 mm^2/s ny90 = 7.25 mm/s^2

rho = 820; % kg/m^3: Average oil density in the operating range

%% Kinematic viscosity of the oil during operation at operating temperature T
% Calculation according to DIN 51563:2011-04
ny40 = 30.8;
T1 = 40;
ny100 = 5.9;
T2 = 100;

W1 = log10(log10(ny40 + 0.8));
W2 = log10(log10(ny100 + 0.8));

m = (W1 - W2) / (log10(T2 + 273.15) - log10(T1 + 273.15));

ny_T = @(T) (10^(10^(m * (log10(T2 + 273.15) - log10(T + 273.15)) + W2))) - 0.8;

ny = ny_T(Temp);


%% Fixed parameters
XL = 0.75; % Lubrication factor
Qe = 2.5 / 1000 / 60; % m^3/s; injection flow rate
v = 0.4; % m/s; injection speed
Ra = 0.5; % µm: Roughness


%%%%%%%%%%%%%%%%%
%% Bearing data
DLager = Gearbox.DLager;
dLager = Gearbox.dLager;
C0Lager = Gearbox.C0Lager;


firstrun = true; % Identifies the first calculation run, no loss values known yet


%% Preallocation
M = zeros(2,1);
Hv = zeros(1,1);
mym = zeros(1,1);
P_VZP = zeros(1,1);
P_VZ0 = zeros(1,1);
P_VQ = zeros(1,1);
P_VI = zeros(1,1);
P_VZ = zeros(1,1);
P_VD = zeros(2,1);
Mrr = zeros(4,1);
Msl = zeros(4,1);
P_VL = zeros(4,1);

Welle = struct('n',{0 0},'P_V',{0 0},'L1Fax',{0 0},'L2Fax',{0 0},'L1Fr',{0 0},'L2Fr',{0 0},'L1Fy',{0 0},'L2Fy',{0 0},'L1Fz',{0 0},'L2Fz',{0 0});
Stufe = struct('P',{0},'P_V',{0},'Fax',{0},'Fu',{0},'Frad',{0},'dw1',{0},'dw2',{0});


%% Three iterative calculation runs (including consideration of
% the losses calculated in the previous run)
for k = 1:3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stage 1 (gear side + wheel side)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    alpha = Gearbox.alpha; % Degrees: Normal pressure angle
    beta = Gearbox.beta(1); % Degrees: Helix angle at the pitch circle
    mn = Gearbox.mn(1); % mm: Normal module
    z1 = Gearbox.z1(1); % Number of teeth on gear 1
    z2 = Gearbox.z2(1); % Number of teeth on gear 2
    b = Gearbox.b(1); % mm: Tooth width
    

    % Calculate gear data from inputs
    % Script for gear data according to Appendix B of the semester project
    % Sources: Niemann, Getriebe allgemein, 2003, p.58 and 276f.
    %          FZG, Formula collection for the Machine Elements module, 2020, p.102ff.
    %          Gruss, Calculation basics for spur gears (available online) 
        alphawt = atand(tand(alpha) / cosd(beta)); % Degrees: Operating pressure angle
        betab = atand(tand(beta) * cosd(alphawt)); % Degrees: Base helix angle
        mt = mn / cosd(beta); % mm: Transverse module
        
        i_Stufe = z2 / z1; % --: Gear ratio; transmission ratio
        
        d1 = mt * z1; % mm: Pitch circle diameter of gear 1
        d2 = mt * z2; % mm: Pitch circle diameter of gear 2
    
        dw1 = d1; % mm: Rolling circle diameter of gear 1
        dw2 = d2; % mm: Rolling circle diameter of gear 2
    
        db1 = d1 * cosd(alphawt); % mm: Base circle diameter of gear 1
        db2 = d2 * cosd(alphawt); % mm: Base circle diameter of gear 2
        
        a = (d1 + d2) / 2; % mm: Center distance
        
        haP = mn; % mm: Addendum
        da1 = d1 + 2 * haP; % mm: Addendum circle diameter of gear 1
        da2 = d2 + 2 * haP; % mm: Addendum circle diameter of gear 2
        
        galpha = 0.5 * (sqrt(da1^2 - db1^2) + sqrt(da2^2 - db2^2)) - a * sind(alphawt); % mm: Path of contact
        galphaa1 = 0.5 * db1 * (((da1 / db1)^2 - 1)^0.5 - tand(alphawt)); % mm: Path of contact for pinion addendum
        eps_alpha = galpha / (mt * pi * cosd(alphawt)); % Profile contact ratio
        eps_1 = galphaa1 / (mt * pi * cosd(alphawt)); % Partial contact ratio for pinion addendum
        eps_2 = eps_alpha - eps_1; % Partial contact ratio for gear addendum
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% P_VZP (Load-dependent gear losses)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Calculate speeds and torques in stage 1
        Welle(2).n = nrad; % Speed of shaft 2
        Welle(1).n = Welle(2).n * i_Stufe; % Speed of shaft 1
        if firstrun % Losses unknown
            M(2) = Mrad; % Torque on shaft 2
        elseif Mrad >= 0 % Driving mode Mrad > 0 -> Torque increases
            M(2) = Mrad + (Stufe(1).P_V + Welle(2).P_V) / (2 * pi * Welle(2).n / 60);
        else % Recuperation mode Mrad < 0 -> Torque decreases
            M(2) = Mrad + Welle(2).P_V / (2 * pi * Welle(2).n / 60);
        end
        
        Stufe(1).P = (Welle(2).n * 2 * pi / 60) * M(2); % Power transmitted in stage 1
        
        % Gear forces
        Fu = 2 * M(2) / (d2 * 1e-3); % Tangential force
        Fax = Fu * tand(beta); % Axial force
        Frad = abs(Fu) * tand(alphawt); % Radial force
        
        % Intermediate variables for P_VZP calculation
        Fbt = abs(Fu / cosd(alphawt)); % Tangential force normal to tooth
        vt = Welle(2).n / 60 * pi * d2 * 1e-3; % Tangential velocity of the gear
        vsigC = 2 * vt * sind(alphawt); % Sliding velocity at the contact point
        etaOil = rho * ny / 1e3; % Dynamic viscosity of the oil
        rhoredC = 0.5 * d1 * sind(alphawt) * i_Stufe(1) / (cosd(betab) * (i_Stufe(1) + 1)); % Equivalent curvature radius at the contact point
        
        % Tooth loss factor Hv according to Ohlendorf
        Hv(1) = pi * (i_Stufe(1) + 1) / (z1 * i_Stufe(1) * cosd(betab)) * (1 - eps_alpha + eps_1^2 + eps_2^2);
        
        % Average gear friction coefficient according to Schlenk
        mym(1) = 0.048 * (Fbt / b / (vsigC * rhoredC))^0.2 * etaOil^-0.05 * XL * Ra^0.25;
        
        % Load-dependent loss power in stage 1
        P_VZP(1) = abs(Stufe(1).P) * mym(1) * Hv(1);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% P_VZ0 (Load-independent gear losses) with injection lubrication
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Squeeze losses according to Mauz
        T_Q = 0; % No squeeze losses if sprayed into the outlet (assumption from WTplus)
        P_VQ(1) = T_Q * Welle(2).n * 2 * pi / 60; % = 0, since T_Q = 0
        
        % Impulse losses according to Ariura
        C1 = 1;
        P_VI(1) = d2 / 2 * 1e-3 * rho * Qe * (vt + v) * C1 * Welle(2).n * 2 * pi / 60;
        
        % Ventilation losses
        % P_VV = 1.37e-9 * vt^1.9 * d1^1.6 * b^0.52 * mn^0.69 * na * 2 * pi / 60
        % Neglected
        
        % Load-independent gear losses
        P_VZ0(1) = P_VQ(1) + P_VI(1);
        
        % Total gear losses
        P_VZ(1) = P_VZP(1) + P_VZ0(1);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Bearing forces on shaft 2 (See force diagram in the appendix of SA)
        a = 80 + b / 2; 
        c = 20 + b / 2; % Additional spacing for differential
        Welle(2).L1Fax = Fax;
        Welle(2).L2Fax = 0;
        Welle(2).L1Fy = (+Frad * c + Fax * dw2 / 2) / (a + c);
        Welle(2).L2Fy = (Frad - Welle(2).L1Fy);
        Welle(2).L2Fz = (Fu * a) / (a + c);
        Welle(2).L1Fz = (Fu - Welle(2).L2Fz);
        Welle(2).L1Fr = sqrt(Welle(2).L1Fy^2 + Welle(2).L1Fz^2);
        Welle(2).L2Fr = sqrt(Welle(2).L2Fy^2 + Welle(2).L2Fz^2);
        
        % Save stage data
        Stufe(1).Fax = Fax;
        Stufe(1).Fu = Fu;
        Stufe(1).Frad = Frad;
        Stufe(1).dw1 = dw1;
        
        
        %% Bearing forces on shaft 1 (See force diagram in the appendix of SA)
        a = 8 + b / 2; 
        c = 8 + b / 2; 
        Welle(1).L1Fax = -Fax;
        Welle(1).L2Fax = 0;
        Welle(1).L1Fy = (Fax * dw1 / 2 - Frad * c) / (a + c);
        Welle(1).L2Fy = -Frad - Welle(1).L1Fy;
        Welle(1).L2Fz = - Fu * a / (a + c);
        Welle(1).L1Fz = - Fu - Welle(1).L2Fz;
        Welle(1).L1Fr = sqrt(Welle(1).L1Fy^2 + Welle(1).L1Fz^2);
        Welle(1).L2Fr = sqrt(Welle(1).L2Fy^2 + Welle(1).L2Fz^2);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% P_VD Seal losses according to ISO 14179-2 (Germany)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Seal at motor input shaft
        P_VD(1) = 7.69e-6 * (Gearbox.d(1)^2 * Welle(1).n);
        
        % Seal at wheel-side output shaft
        P_VD(2) = 2 * 7.69e-6 * (Gearbox.d(2)^2 * Welle(2).n);     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VL Bearing losses SKF2020 (Hardcoded for 160/161 deep groove ball bearings)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Formulas according to SKF, The SKF model for calculation of the frictional moment
    % Online: https://cdn.skfmediahub.skf.com/api/public/0901d196809bc183/pdf_preview_medium/0901d196809bc183_pdf_preview_medium.pdf
    
    % Bearing rotational speeds = Shaft rotational speeds
    nLager = [Welle(1).n Welle(2).n];
    
    for i = 1:2 % Iteration over shafts
        for j = 1:2 % Iteration over bearing positions per shaft
            index = 2 * (i - 1) + j;
    
            % Read specific bearing values from the design
            DLag = DLager(i); % Outer diameter of the bearing
            dLag = dLager(i); % Inner diameter of the bearing
            C0 = C0Lager(i); % Static load rating of the bearing
            nLag = nLager(i); % Bearing rotational speed
            
            dm = (DLag + dLag) / 2; % Mean bearing diameter
    
            % Assign acting forces
            if j == 1 % Bearing 1
                Fr = Welle(i).L1Fr; Fa = abs(Welle(i).L1Fax);
            else % Bearing 2
                Fr = Welle(i).L2Fr; Fa = abs(Welle(i).L2Fax);
            end
            
            % Coefficients for 160/161 deep groove ball bearings
            S1 = 4.63e-3; S2 = 4.25;
            R1 = 4.3e-7; R2 = 1.7;
            
            % Lubrication film thickness factor psi_ish
            psiish = 1 / (1 + 1.84e-9 * (nLag * dm)^1.28 * ny^0.64);
    
            % Lubricant displacement factor psi_rs
            psirs = 1 / exp(6e-8 * ny * nLag * (DLag + dLag) * sqrt(3.1 / (2 * (DLag - dLag))));
    
            % Rolling friction base value G_rr, Sliding friction base value G_sl
            if Fa == 0 % Without axial load
                Grr = R1 * dm^1.96 * Fr^0.54;
                Gsl = S1 * dm^(-0.26) * Fr^(5 / 3);
            else % With axial load
                alphaf = 24.6 * (Fa / C0)^0.24;
                Grr = R1 * dm^1.96 * (Fr + R2 / sind(alphaf) * Fa)^0.54;
                Gsl = S1 * dm^-0.145 * (Fr^5 + S2 * dm^1.5 / sind(alphaf) * Fa^4)^(1 / 3);
            end
    
            % Rolling friction moment M_rr
            Mrr(index) = psiish * psirs * Grr * (ny * nLag)^0.6;
            
            % Weighting factor for sliding friction coefficient psi_bl
            psibl = 1 / exp(2.6e-8 * (nLag * ny)^1.4 * dm);
    
            % Sliding friction coefficient my_sl
            mysl = psibl * 0.12 + (1 - psibl) * 0.04; % SKF online calculator uses 0.05 instead of 0.04
    
            % Sliding friction moment M_sl
            Msl(index) = Gsl * mysl;
            
            % Bearing power loss
            P_VL(index) = (Mrr(index) + Msl(index)) * 1e-3 * nLag * 2 * pi / 60;
        end
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Assignment of losses to the corresponding shaft or stage
    % Necessary to account for them correctly in the next iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Stufe(1).P_V = P_VZ(1); % Assign gear losses to stage 1
    Welle(1).P_V = P_VD(1) + sum(P_VL(1:2)); % Assign seal and bearing losses to shaft 1
    Welle(2).P_V = P_VD(2) + sum(P_VL(3:4)); % Assign seal and bearing losses to shaft 2
    
    firstrun = false; % In subsequent iterations, losses are known
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Total loss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_V = sum(P_VD) + sum(P_VZ) + sum(P_VL); % Calculate total gearbox loss

% Uncomment to display the loss components in the Command Window:
% disp(['Total gear loss power: ', sprintf('%0.2f', sum(P_VZ)), ' W'])
% disp(['Total seal loss power: ', sprintf('%0.2f', sum(P_VD)), ' W'])
% disp(['Total bearing loss power: ', sprintf('%0.2f', sum(P_VL)), ' W'])
% disp(['Total gearbox loss power: ', sprintf('%0.2f', sum(P_VD) + sum(P_VZ) + sum(P_VL)), ' W'])

end