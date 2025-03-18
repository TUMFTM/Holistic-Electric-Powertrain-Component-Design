%Skript für Verzahnungsdaten
% Quellen: Niemann, Getriebe allgemein, 2003, S.58 und 276f.
%          FZG, Formelsammlung zum Modul Maschinenelemente, 2020, S.102ff.
%          Gruss, Berechnungsgrundlagen Stirnradverzahnungen (online verfügbar) 

    alphawt = atand(tand(alpha)/cosd(beta)); %Grad: Betriebseingriffswinkel
    betab = atand(tand(beta)*cosd(alphawt)); %Grad: Schrägungswinkel Grundkreis
    mt = mn/cosd(beta); %mm: Stirnmodul
    
    u = z2/z1; %--: Zähneverhältnis; Übersetzung
    
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
