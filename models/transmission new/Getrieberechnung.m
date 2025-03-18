function [P_V] = Getrieberechnung(Gearbox,Mrad,nrad,Temp)
% Funktion zur Verlustberechnung eines 2-stufigen Strinradgetriebes mit
% achsparallelem An/Abtrieb, ebene Anordnung

%%#codegen

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Öl- und Schmierungsdaten
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%%% Öldaten SAE 70W-75W (API GL4) VW ID3 Getriebeöl nach Liqui Moly Datenblatt
% ny40 = 30.8 mm^2/s  und     ny100 = 5.9 mm^2/s    rho15 = 845 kg/m^3
% ny80 = 9.12 mm^2/s  ny90 = 7.25 mm/s^2

rho = 820; %kg/m^3: mittlere Dichte des Öls im Betriebsbereich

%% kinematische Viskosität des Öls im Betrieb bei Betriebstemperatur T
%Berechnung nach DIN 51563:2011-04
ny40 = 30.8;
T1 = 40;
ny100 = 5.9;
T2 = 100;

W1 = log10(log10(ny40+0.8));
W2 = log10(log10(ny100+0.8));

m = (W1-W2)/(log10(T2+273.15)-log10(T1+273.15));

ny_T = @(T) (10^(10^(m*(log10(T2+273.15)-log10(T+273.15))+W2)))-0.8;

ny = ny_T(Temp);


%% festgelegte Parameter
XL = 0.75; % Schmierstofffaktor
Qe = 2.5/1000/60; %m^3/s; Einspritzvolumenstrom
v = 0.4; %m/s; Einspritzgeschwindigkeit
Ra = 0.5; %ym: Rauigkeit


%%%%%%%%%%%%%%%%%
%% Lagerdaten
DLager = Gearbox.DLager;
dLager = Gearbox.dLager;
C0Lager = Gearbox.C0Lager;


firstrun = true; % Identifiziert den ersten Berechnungsdurchgang, noch keine Verlustwerte bekannt


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


%% Dreifacher, iterativer Berechnungsdurchlauf (inkl. Berücksichtigung der
% im vorherigen Durchgang berechneten Verlustwerte)
for k = 1:3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stufe 2 (radseitig)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha = Gearbox.alpha; %Grad: Normaleingriffswinkel
    beta = Gearbox.beta(2); %Grad: Schrägungswinkel Teilkreis
    mn = Gearbox.mn(2); %mm: Normalmodul
    z1 = Gearbox.z1(2); %Zähnezahl Zahnrad 1
    z2 = Gearbox.z2(2); %Zähnezahl Zahnrad 2
    b = Gearbox.b(2); %mm: Zahnbreite
    

    % Verzahnungsdaten aus Eingaben berechnen
    %Skript für Verzahnungsdaten nach Anhang B der Semesterarbeit
    % Quellen: Niemann, Getriebe allgemein, 2003, S.58 und 276f.
    %          FZG, Formelsammlung zum Modul Maschinenelemente, 2020, S.102ff.
    %          Gruss, Berechnungsgrundlagen Stirnradverzahnungen (online verfügbar) 
        alphawt = atand(tand(alpha)/cosd(beta)); %Grad: Betriebseingriffswinkel
        betab = atand(tand(beta)*cosd(alphawt)); %Grad: Schrägungswinkel Grundkreis
        mt = mn/cosd(beta); %mm: Stirnmodul
        
        i_Stufe = z2/z1; %--: Zähneverhältnis; Übersetzung
        
        d1 = mt*z1; %mm: Teilkreisdurchmesser Zahnrad 1
        d2 = mt*z2; %mm: Teilkreisdurchmesser Zahnrad 2
    
        dw1 = d1; %mm: Wälzkreisdurchmesser Z1
        dw2 = d2; %mm: Wälzkreisdurchmesser Z2
    
        db1 = d1 * cosd(alphawt); %mm: Grundkreisdurchmesser Z1
        db2 = d2 * cosd(alphawt); %mm: Grundkreisdurchmesser Z2
        
        a = (d1 +d2)/2; %mm: Achsabstand
        
        haP = mn; %mm: Kopfhöhe
        da1 = d1 + 2*haP; %mm: Kopfkreisdurchmesser Z1
        da2 = d2 + 2*haP; %mm: Kopfkreisdurchmesser Z2
        
        galpha = 0.5* (sqrt(da1^2 - db1^2) + sqrt(da2^2 - db2^2)) - a*sind(alphawt); %mm: Eingriffsstrecke
        galphaa1 = 0.5*db1 * (((da1/db1)^2 -1)^0.5 - tand(alphawt)); %mm: Eingriffsstrecke Ritzelkopf
        eps_alpha = galpha / (mt*pi*cosd(alphawt)); %Profilüberdeckung
        eps_1 = galphaa1 / (mt*pi*cosd(alphawt)); %Teilüberdeckung Ritzelkopf
        eps_2 = eps_alpha - eps_1; %Teilüberdeckung Radkopf
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZP (lastabhängige Verzahnungsverluste)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Drehzahlen und Drehmomente in Stufe 2 berechnen    
    Welle(3).n = nrad;
    Welle(2).n = Welle(3).n*i_Stufe;
    if firstrun %Verluste unbekannt
        M(3) = Mrad; 
    elseif Mrad >= 0 %Antriebsfall Mrad>0 -> Moment erhöht
        M(3) = Mrad + (Stufe(2).P_V+Welle(3).P_V)/(2*pi*Welle(3).n/60);
    else %Rekuperationsfall Mrad<0 -> Moment reduziert
        M(3) = Mrad + Welle(3).P_V/(2*pi*Welle(3).n/60);
    end
    
    Stufe(2).P = (Welle(3).n*2*pi/60)*M(3); %übertragene Leistung in Stufe 2
    
    
    % Verzahnungskräfte
    Fu = 2*M(3)/(d2*1e-3); %Umfangskraft
    Fax = Fu * tand(beta); %Axialkraft
    Frad = abs(Fu) * tand(alphawt); %Radialkraft
    
    % Zwischengrößen für P_VZP-Berechnung
    Fbt = abs(Fu/cosd(alphawt)); %Umfangskraft normal auf Zahn
    vt = Welle(3).n/60 * pi * d2*1e-3; %Umfangsgeschwindigkeit Zahnrad
    vsigC = 2*vt*sind(alphawt); %Summengeschwindigkeit Wälzpunkt
    etaOil = ny*rho/1e3; %dynamische Viskosität des Öls
    rhoredC = 0.5*d1*sind(alphawt)*i_Stufe(1)/(cosd(betab)*(i_Stufe(1)+1)); %Ersatzkrümmungsradius Wälzpunkt
    
    %Zahnverlustfaktor Hv nach Ohlendorf
    Hv(2) = pi*(i_Stufe(1)+1)/(z1*i_Stufe(1)*cosd(betab))*(1-eps_alpha + eps_1^2+eps_2^2);
    
    %Mittlere Verzahnungsreibungszahl nach Schlenk
    mym(2) =0.048*(Fbt/b/(vsigC*rhoredC))^0.2 * etaOil^-0.05*XL*Ra^0.25;
    
    %lastabhängige Verlustleistung in Stufe 2
    P_VZP(2) = abs(Stufe(2).P) * mym(2) * Hv(2);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZ0 (lastunabhängige Verzahnungsverluste) mit Einspritzschmierung
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Quetschverluste nach Mauz
    T_Q = 0; %Keine Quetschverluste, wenn in Auslauf gespritzt wird (Annahme aus WTplus)
    P_VQ(2) = T_Q * Welle(2).n *2*pi/60; %= 0, da T_Q = 0
    
    %Impulsverluste nach Ariura
    C1 = 1; 
    P_VI(2) = d2/2 *1e-3 * rho * Qe * (vt+v) * C1 * Welle(3).n *2*pi/60;
    
    %Ventilationsverluste
    %P_VV = 1.37e-9 * vt^1.9 * d1^1.6 * b^0.52 * mn^0.69 * na *2*pi/60
    % vernachlässigt
    
    %lastunabhängige Verzahnungsverluste
    P_VZ0(2) = P_VQ(2) + P_VI(2);

    %Summe der Verzahnungsverluste
    P_VZ(2) = P_VZP(2) + P_VZ0(2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Lagerkräfte Welle 3 (Siehe Kraftplan im Anhang der SA)
    a = 80+b/2; b1 = 0; c = 20+b/2; %Extraabstand für Differential
    Welle(3).L1Fax = Fax;
    Welle(3).L2Fax = 0;
    Welle(3).L1Fy = -(Frad*c + Fax*dw2/2)/(a+b1+c);
    Welle(3).L2Fy = -(Frad + Welle(3).L1Fy);
    Welle(3).L2Fz = (Fu*a)/(a+b1+c);
    Welle(3).L1Fz = (Fu - Welle(3).L2Fz);
    Welle(3).L1Fr = sqrt(Welle(3).L1Fy^2 + Welle(3).L1Fz^2);
    Welle(3).L2Fr = sqrt(Welle(3).L2Fy^2 + Welle(3).L2Fz^2);
    
    %Abspeichern der Stufendaten
    Stufe(2).Fax = Fax;
    Stufe(2).Fu = Fu;
    Stufe(2).Frad = Frad;
    Stufe(2).dw1 = dw1;
    
    % Berechnung des Drehmoments an Welle 2
    if firstrun %Verluste unbekannt
        M(2) = M(3)/i_Stufe;
    elseif Mrad >= 0 %Antriebsfall, Verluste erhöhen Leistung
        M(2) = M(3)/i_Stufe + (Stufe(1).P_V+Welle(2).P_V)/(2*pi*Welle(2).n/60);
    else %Rekuperationsfall, Verluste redzieren negative Leistung
        M(2) = M(3)/i_Stufe + (Stufe(2).P_V+Welle(2).P_V)/(2*pi*Welle(2).n/60);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stufe 1 (motorseitig)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = Gearbox.beta(1); %Grad: Schrägungswinkel Teilkreis
    mn = Gearbox.mn(1); %mm: Normalmodul
    z1 = Gearbox.z1(1); %Zähnezahl Zahnrad 1
    z2 = Gearbox.z2(1); %Zähnezahl Zahnrad 2
    b = Gearbox.b(1); %mm: Zahnbreite
    
    %Verzahnungsdaten aus Eingaben berechnen
    %Skript für Verzahnungsdaten nach Anhang B der Semesterarbeit
    % Quellen: Niemann, Getriebe allgemein, 2003, S.58 und 276f.
    %          FZG, Formelsammlung zum Modul Maschinenelemente, 2020, S.102ff.
    %          Gruss, Berechnungsgrundlagen Stirnradverzahnungen (online verfügbar) 
        alphawt = atand(tand(alpha)/cosd(beta)); %Grad: Betriebseingriffswinkel
        betab = atand(tand(beta)*cosd(alphawt)); %Grad: Schrägungswinkel Grundkreis
        mt = mn/cosd(beta); %mm: Stirnmodul
        
        i_Stufe = z2/z1; %--: Zähneverhältnis; Übersetzung
        
        d1 = mt*z1; %mm: Teilkreisdurchmesser Zahnrad 1
        d2 = mt*z2; %mm: Teilkreisdurchmesser Zahnrad 2
    
        dw1 = d1; %mm: Wälzkreisdurchmesser Z1
        dw2 = d2; %mm: Wälzkreisdurchmesser Z2
    
        db1 = d1 * cosd(alphawt); %mm: Grundkreisdurchmesser Z1
        db2 = d2 * cosd(alphawt); %mm: Grundkreisdurchmesser Z2
        
        a = (d1 +d2)/2; %mm: Achsabstand
        
        haP = mn; %mm: Kopfhöhe
        da1 = d1 + 2*haP; %mm: Kopfkreisdurchmesser Z1
        da2 = d2 + 2*haP; %mm: Kopfkreisdurchmesser Z2
        
        galpha = 0.5* (sqrt(da1^2 - db1^2) + sqrt(da2^2 - db2^2)) - a*sind(alphawt); %mm: Eingriffsstrecke
        galphaa1 = 0.5*db1 * (((da1/db1)^2 -1)^0.5 - tand(alphawt)); %mm: Eingriffsstrecke Ritzelkopf
        eps_alpha = galpha / (mt*pi*cosd(alphawt)); %Profilüberdeckung
        eps_1 = galphaa1 / (mt*pi*cosd(alphawt)); %Teilüberdeckung Ritzelkopf
        eps_2 = eps_alpha - eps_1; %Teilüberdeckung Radkopf
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZP (lastabhängige Verzahnungsverluste) in Stufe 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Drehzahl Welle 1 und übertragende Leistung Stufe 1
    Welle(1).n = Welle(2).n*i_Stufe;
    Stufe(1).P = (Welle(2).n*2*pi/60)*M(2);
    
    %Verzahnungskräfte
    Fu = 2*M(2)/(d2*1e-3); %Umfangskraft
    Fax = -Fu * tand(beta); %Axialkraft
    Frad = abs(Fu)*tand(alphawt); %Radialkraft
    
    %Zwischengrößen für P_VZP-Berechnung
    Fbt = abs(Fu/cosd(alphawt)); %Umfangskraft normal auf Zahn
    vt = Welle(2).n/60 * pi * d2*1e-3; %Umfangsgeschwindigkeit Zahnrad
    vsigC = 2*vt*sind(alphawt); %Summengeschwindigkeit Wälzpunkt
    etaOil = rho*ny/1e3; %dynamische Viskosität Öl
    rhoredC = 0.5*d1*sind(alphawt)*i_Stufe(1)/(cosd(betab)*(i_Stufe(1)+1)); %Ersatzkrümmungsradius Wälzpunkt
    
    %Zahnverlustfaktor Hv nach Ohlendorf
    Hv(1) = pi*(i_Stufe(1)+1)/(z1*i_Stufe(1)*cosd(betab))*(1-eps_alpha + eps_1^2+eps_2^2);
    
    %Mittlere Verzahnungsreibungszahl nach Schlenk
    mym(1) =0.048*(Fbt/b/(vsigC*rhoredC))^0.2 * etaOil^-0.05*XL*Ra^0.25;
    
    %lastabhängige Verlustleistung in Stufe 1
    P_VZP(1) = abs(Stufe(1).P) * mym(1) * Hv(1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VZ0 (lastunabhängige Verzahnungsverluste) mit Einspritzschmierung
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %Quetschverluste nach Mauz
    T_Q = 0; %Keine Quetschverluste, wenn in Auslauf gespritzt wird (WTplus Annahme)
    P_VQ(1) = T_Q * Welle(1).n *2*pi/60;
    
    % Impulsverluste nach Ariura
    C1 = 0.85; 
    P_VI(1) = d2/2*1e-3 * rho * Qe * (vt+v) * C1 * Welle(2).n *2*pi/60;

    % Ventilationsverluste P_VV vernachlässigt
    
    % lastunabhängige Verzahnungsverluste Stufe 1:
    P_VZ0(1) = P_VQ(1) + P_VI(1);

    % Verzahnungsverluste Stufe 1 gesamt:
    P_VZ(1) = P_VZP(1) + P_VZ0(1);
    
    % Abspeichern der Stufendaten
    Stufe(1).Fax = Fax;
    Stufe(1).Fu = Fu;
    Stufe(1).Frad = Frad;
    Stufe(1).dw2 = dw2;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Lagerkräfte Welle 2 (siehe Kräfteplan im Anhang der SA)
    a = 8+Gearbox.b(1)/2; b1 = Gearbox.b(1)/2+5+Gearbox.b(2)/2; c = Gearbox.b(2)/2+8;
    Welle(2).L1Fax = - Stufe(1).Fax - Stufe(2).Fax;
    Welle(2).L2Fax = 0;
    Welle(2).L1Fy = (Stufe(1).Fax*Stufe(1).dw2/2 - Stufe(1).Frad*(b1+c) + Stufe(2).Frad*c - Stufe(2).Fax*Stufe(2).dw1/2)/(a+b1+c);
    Welle(2).L2Fy = -(Stufe(1).Frad - Stufe(2).Frad + Welle(2).L1Fy);
    Welle(2).L2Fz = -(Stufe(1).Fu*a + Stufe(2).Fu*(a+b1)) / (a+b1+c);
    Welle(2).L1Fz = -(Stufe(1).Fu + Stufe(2).Fu + Welle(2).L2Fz);
    Welle(2).L1Fr = sqrt(Welle(2).L1Fy^2 + Welle(2).L1Fz^2);
    Welle(2).L2Fr = sqrt(Welle(2).L2Fy^2 + Welle(2).L2Fz^2);
    
    
    %% Lagerkräfte Welle 1 (siehe Kräfteplan im Anhang der SA)
    a = 8+b/2; b1 = 0; c = b/2+8; 
    Welle(1).L1Fax = Fax;
    Welle(1).L2Fax = 0;
    Welle(1).L1Fy = (Fax*dw1/2 + Frad*c)/(a+b1+c);
    Welle(1).L2Fy = Frad - Welle(1).L1Fy;
    Welle(1).L2Fz = Fu*a/(a+b1+c);
    Welle(1).L1Fz = Fu- Welle(1).L2Fz;
    Welle(1).L1Fr = sqrt(Welle(1).L1Fy^2 + Welle(1).L1Fz^2);
    Welle(1).L2Fr = sqrt(Welle(1).L2Fy^2 + Welle(1).L2Fz^2);
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VD Dichtungsverluste nach ISO 14179-2 (Deutschland)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Dichtung Eingangswelle Motor
    P_VD(1) = 7.69e-6 * (Gearbox.d(1)^2 * Welle(1).n);
    
    %Dichtungen (2x) Ausgangswelle Rad
    P_VD(2) = 2* 7.69e-6 * (Gearbox.d(3)^2 * Welle(3).n);

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% P_VL Lagerverluste SKF2020 (Hardcoded für 160/161er Rillenkugellager)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Formeln nach SKF, The SKF model for calculation the frictional moment
    % Online: https://cdn.skfmediahub.skf.com/api/public/0901d196809bc183/pdf_preview_medium/0901d196809bc183_pdf_preview_medium.pdf

    %Lagerdrehzahlen = Wellendrehzahlen
    nLager = [Welle(1).n Welle(2).n Welle(3).n];
    
    for i  = 1:3 %Iteration über Wellen
        for j = 1:2 %Iteration über Lagerstellen je Welle
            index = 2*(i-1)+j;

            %spez. Lagerwerte aus Design einlesen
            DLag = DLager(i); %Außendurchmesser Lager
            dLag = dLager(i); %Innendurchmesser Lager
            C0 = C0Lager(i); %stat. Tragzahl Lager
            nLag = nLager(i); %Lagerdrehzahl
            
            dm = (DLag + dLag)/2; %mittlerer Lagerdurchmesser
            
            %Zuordnen der wirkenden Kräfte
            if j == 1 %Lager 1
                Fr = Welle(i).L1Fr; Fa = abs(Welle(i).L1Fax);
            else %Lager 2
                Fr = Welle(i).L2Fr; Fa = abs(Welle(i).L2Fax);
            end
            
            %Koeffs für Rillenkugellager 160/161er Reihe
            S1 = 4.63e-3; S2 = 4.25;
            R1 = 4.3e-7; R2 = 1.7;
            
            %Schmierfilmdickenfaktor psi_ish
            psiish = 1/(1 + 1.84e-9 * (nLag*dm)^1.28 * ny^0.64);

            %Schmierstoffverdrängungsfaktor psi_rs
            psirs = 1/exp(6e-8 * ny * nLag * (DLag+dLag) * sqrt(3.1/(2*(DLag - dLag))));

            %Rollreibungsgrundwert G_rr, Gleitreibungsgrundwert G_sl
            if Fa == 0 %ohne Axiallast
                Grr = R1 * dm^1.96 * Fr^0.54;
                Gsl = S1 * dm^(-0.26)*Fr^(5/3);
            else %bei Axialbelastung
                alphaf = 24.6*(Fa/C0)^0.24;
                Grr = R1 * dm^1.96 * (Fr + R2/sind(alphaf)*Fa)^0.54;
                Gsl = S1 * dm^-0.145 * (Fr^5 + S2*dm^1.5/sind(alphaf)*Fa^4)^(1/3);
            end

            %Rollreibungsmoment M_rr
            Mrr(index) = psiish * psirs * Grr * (ny*nLag)^0.6;
            
            %Gewichtungsfaktor Gleitreibungszahl psi_bl
            psibl = 1/exp(2.6e-8 * (nLag*ny)^1.4 * dm);

            %Gleitreibungszahl my_sl
            mysl = psibl*0.12 + (1-psibl)*0.04; %SKF Online-Rechner mit 0.05 statt 0.04

            %Gleitreibungsmoment M_sl
            Msl(index) = Gsl * mysl;
            
            %Lagerverlustleistung
            P_VL(index) = (Mrr(index) + Msl(index)) *1e-3 * nLag*2*pi/60;
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Zuordnung der Verluste zu entsprechender Welle oder Stufe
    % notwendig, um sie im nächsten Durchgang an der richtigen Stelle zu berücksichtigen
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Stufe(1).P_V = P_VZ(1);
    Stufe(2).P_V = P_VZ(2);
    Welle(1).P_V = P_VD(1) + sum(P_VL(1:2));
    Welle(2).P_V = sum(P_VL(3:4));
    Welle(3).P_V = P_VD(2) + sum(P_VL(5:6));
    
    firstrun = false; % in nachfolgenden Durchgängen sind Verluste bekannt
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gesamtverlust
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_V = sum(P_VD)+sum(P_VZ)+sum(P_VL);

% Zur Ausgabe der Verlustanteile im Command Window aktivieren:
% disp(['Verzahnungsverlustleistung Gesamt: ',sprintf('%0.2f',sum(P_VZ)),' W'])
% disp(['Dichtungsverlustleistung Gesamt: ',sprintf('%0.2f',sum(P_VD)),' W'])
% disp(['Lagerverlustleistung Gesamt: ',sprintf('%0.2f',sum(P_VL)),' W'])
% disp(['Getriebeverlustleistung Gesamt: ',sprintf('%0.2f',P_V),' W'])

end