function [Gearbox] = Getriebedesign_koaxial(M,iG)
% Funktion zum Entwurf eines koaxialen Getriebes (zweistufig)
% Input: ertragbares Drehmoment M, 
%        Gesamtübersetzung iG, 
Mmax = M;
iges = iG;


%% Berechnung

simga_H_lim = 1500; %N/mm^2, Dauer/Zeitfestigkeitswert Naunheimer, Fahrzeuggetriebe, 2019, Tab. 7.1

alpha = 20; %Grad; Normaleingriffswinkel; 20° = Normwert
beta = [23.4 18.4]; %Grad; Schrägungswinkel; übernommen von Gao
AnzahlStufen = 2; %muss zweistufig sein
    
%Aufteilung der Übersetzung auf beide Stufen, sodass a(1) == a(2) für Koaxialität
%Formel für Achsabstand nach Naunheimer, Fahrzeuggetriebe, 2019, S.288
syms i1; %symbolische Variable für Übersetzung Stufe 1
assume(i1,'real');
assumeAlso(i1>0);
M(1) = Mmax;
eq = 0.299^3*M(1)*(i1+1)^4/i1 == 0.255^3*M(1)*i1*(iges/i1+1)^4/(iges/i1); %a(1) == a(2)
sol = solve(eq,i1);
i= [double(sol)];
i = [i iges/i]; %Übersetzungen Stufe 1 und 2

% Erfahrungswerte Ritzelzähnezahl aus Niemann Getriebe allgemein, 2003, S.269, einsatzgehärtet
i_inter = [1 2 4 8];
z_inter1 = [32 29 25 22]; %hohe Drehzahl
z_inter2 = [26.5 24 20.5 18]; %mittlere Drehzahl

% Interpolation der Zähnezahl für Stufe 1 und 2 basierend auf Erfahrungswerten
z(1) = round(interp1(i_inter,z_inter1,i(1),'linear'),0); %hohe Drehzahlen
z(2) = round(interp1(i_inter,z_inter2,i(2),'linear'),0); %mittlere Drehzahlen

% Berechnung der Radzähnezahl(en)
z2 = round(z.*i);


%% Wellendurchmesser
% Wellendurchmesser für Dichtringgröße und Lagerauswahl notwendig

%Vergleichsmoment, da Biege- und Torsionsbelastung: Faktor 1.2 ... 2.5
%Quelle Niemann Winter Maschinenelemente 1, 2019, Formel 17.10
%Mv = 1.5*Mmax; -> unten eingesetzt

%Stahlannahme: 16MnCr5 sigma_bW = 520 N/mm^2 aus Naunheimer, Fahrzeuggetriebe, 2019, S.478
%Sicherheitsfaktor 1,5   Größeneinfluss 0.85   Oberflächeneinfluss 0.8

sigma_b_zul = 520*0.85*0.8/1.5; %zulässige Biegespannung nach Naunheimer, Fahrzeuggetriebe, 2019, S.478

for j=1:AnzahlStufen+1
    Mv = 1.5*Mmax*prod(i(1:j-1)); %Vergleichsmoment
    if j > AnzahlStufen %wenn Abtriebswele
        Mv = Mv/2; %Auteilung des Moments auf links + rechts
    end
    d(j) = 2.17 * nthroot(Mv*1e3/sigma_b_zul,3); %Wellendurchmesser nach Naunheimer, Fahrzeuggetriebe, 2019, S.476
end


%Durchmesser Zwischenwelle nach Niemann, Maschinenelemente 1, 2019, S.494
%zulässige Schubspannung 190 N/mm^2 basierend auf Taycan-Welle ermittelt
dZwischenwelle = 1.72 * nthroot(0.5*Mmax*1e3*iges/190,3);


%Außendurchmesser Hohlwelle nach Wittel, Roloff/Matek Maschinenelemente, 2021, S. 388
syms d1a; %symbolische Variable für Hohlwellenaußendurchmesser
assume(d1a,'real');
assumeAlso(d1a>0);
eq = 32 * 1.5*Mmax*1e3 / (pi * sigma_b_zul) == (d1a^4 - (dZwischenwelle+1)^4)/d1a;
sol = solve(eq,d1a);
d(1) = double(sol); %Außendurchmesser Hohlwelle


%% Achsabstand a, Zahnbreite b und Normalmodul mn

M(1) = Mmax; %Drehmoment Stufe 1
M(2) = M(1)*i(1); %Drehmoment Stufe 2

for j = 1:AnzahlStufen
    %Achsabstand nach Naunheimer, Fahrzeuggetriebe, 2019, S.288 
    if j==1 %Stufe 1
        a(j) = 0.299 * nthroot((M(j)*1e3*(i(j)+1)^4)/i(j),3);
    else %Stufe 2
        a(j) = 0.255 * nthroot((M(j)*1e3*(i(j)+1)^4)/i(j),3);
    end
    
    % Wälzkreisdurchmesser nach Niemann Getriebe allgemein, 2003, S.276
    d1(j) = 2*a(j)/(1+z2(j)/z(j));
    d2(j) = 2*a(j)/(1+z(j)/z2(j));
    
    %Zahnbreite nach nach Naunheimer, Fahrzeuggetriebe, 2019, S.289
    b(j) = 4278e2 * 0.65 * (M(j)*1e3*(i(j)+1))/(d1(j)^2*i(j)*simga_H_lim^2);

    %Normalmodul nach Wittel, Roloff/Matek Maschinenelemente, 2021, S. 813
    mn(j) = (2*a(j)*cosd(beta(j)))/(z(j)+z2(j));

    %disp(['Pitch circle diameter stage ',int2str(j),': ',sprintf('%0.2f  %0.2f',d1(j),d2(j)),' mm'])
    %disp(['Tooth width stage ',int2str(j),': ',sprintf('%0.2f',b(j)),' mm'])
end


%% Warnung bei geometrischen Problemen mit Hohlwelle und Verzahnung an Eingangswelle

if d1(1) - 2*mn(1) < d(1) %Check Fußkreisdurchmesser < Außendurchmesser Hohlwelle
    disp("- -- -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    disp("Warning: Geometric problem: Hollow shaft to big for calculated pinion in stage 1")
    disp(['Outer diameter hollow shaft: ',sprintf('%0.2f',d(1)),' mm'])
    disp(['Base circle diameter pinion: ',sprintf('%0.2f',d1(1)-2*mn(1)),' mm'])
    disp("- -- -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
end


%% Lagerauswahl

%SKF Lagerkatalog Reihe 160: 16002 bis 16020
di = [15 17 20:5:100]; %Innendurchmesser
da = [32 35 42 47 55 62 68 75 80 90 95 100 110 115 125 130 140 145 150]; %Außendurchmesser
c_stat = [2.85 3.25 4.05 4.75 7.35 8.15 10.2 10.8 11.4 14 15 19.6 25 27 31.5 33.5 39 41.5 44]*1e3; %stat. Tragzahl

for j = 1:AnzahlStufen+1 %Iteration über Wellen, beide Lager je Welle identisch
    %nächstgrößeres Lager auswählen
    Lager_index = sum(d(j) > di)+1;
    % Werte aus Lagerkatalog zuweisen
    dLager(j) = di(Lager_index);
    DLager(j) = da(Lager_index);
    C0Lager(j) = c_stat(Lager_index);
end


%% Abspeichern des Designs im Gearbox-Struct
Gearbox.Stufen = AnzahlStufen; %Stufenanzahl
Gearbox.iges = iges; %Gesamtübersetzung
Gearbox.alpha = alpha; %Normaleingriffswinkel
Gearbox.beta = beta; %Schrägungswinkel
Gearbox.z1 = z; %Zähnezahl(en) Ritzel
Gearbox.z2 = z2; %Zähnezahl(en) Rad/Räder
Gearbox.d1 = d1; %Wälzkreisdurchmesser Ritzel
Gearbox.d2 = d2; %Wälzkreisdurchmesser Rad/Räder
Gearbox.a = a; %Achsabstand/-stände
Gearbox.b = b; %Zahnbreite(n)
Gearbox.mn = mn; %Normalmodul(n)
Gearbox.d = d; %Wellendurchmesser
Gearbox.dLager = dLager; %Lagerinnendurchmesser
Gearbox.DLager = DLager; %Lageraußendurchmesser
Gearbox.C0Lager = C0Lager; %stat. Tragzahlen Lager
end

