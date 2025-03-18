function [Gearbox] = Getriebedesign(M,iG,AnzahlStufen)
% Funktion zum Entwurf eines achsparallelen Getriebes, ein- oder zweistufig
% Input: ertragbares Drehmoment M, 
%        Gesamtübersetzung iG, 
%        Anzahl der Stufen AnzahlStufen (1 oder 2)

Mmax = M;
iges = iG;

%% Refernezwerte aus A2Mac1 für VW ID.3 und Tesla Model 3
% ID3 Getriebe aus A2Mac1 (i1 = 2.95)  (i2 = 3.9)
% Stufe 1: Ritzel: z=23   dk=50mm   b=35mm
%          Rad: z=68   dk=135mm   b=32mm
%
% Stufe 2: Ritzel: z=20   dk=63mm   b=49mm
%          Rad: z=78   dk=217mm   b=45mm

% Tesla M3 Getriebe aus A2Mac1  (i1 = 2.61)  (i2 = 3.46)
% Stufe 1: Ritzel: z=31   dk=58mm   b=37mm
%          Rad: z=81   dk=142mm   b=37mm
%
% Stufe 2: Ritzel: z=24   dk=67mm   b=46mm
%          Rad: z=83   dk=212mm   b=46mm


%% Berechnung

simga_H_lim = 1500; %N/mm^2, Dauer/Zeitfestigkeitswert Naunheimer, Fahrzeuggetriebe, 2019, Tab. 7.1

alpha = 20; %Grad; Normaleingriffswinkel; 20° = Normwert
beta = [23.4 18.4]; %Grad; Schrägungswinkel; übernommen von Gao

% Literaturempfehlungen Stufenanzahl
% Einstufig bis i=5, ansonsten Zweistufig nach Sauer, Konstruktionselemente, 2018, S.480 
% Einstufig bis i=6, ansonsten Zweistufig nach Niemann, Getriebe allgemein, 2003, S.261

if AnzahlStufen > 1 %Zweistufig, erfordert Aufteilung der Übersetzung auf Stufe 1 und 2
    
    %Aufteilung der Übersetzung auf Stufe 1 und 2 nach Niemann Getriebe allgemein, 2003, S.261
    % i(1) = 0.8*iges^(2/3);
    % i(2) = iges/i(1);
    
    %Aufteilung der Übersetzung Stufe 1 und 2 nach Sauer,
    %Konstruktionselemente des Maschinenbaus 2, Seite 480 für 5<iges<15
    %basierend auf Römhild Iris 1993
    i(1) = 0.7332*iges^0.6438;
    i(2) = iges/i(1);
    
    %Beachtung des integrierten Differentials in letzte Stufe fordert ggf.
    %eine größere Übersetzung in Stufe 2, um Platz für Diff zu schaffen
    %Folgende Zeilen können als Overwrite aktiviert werden:
    % i(2) = [3.9]; %Overwrite ID3 [2.95 3.9]
    % i(2) = [3.46]; %Overwrite Tesla M3 [2.61 3.46]
    i(1) = iges/i(2);


    % Erfahrungswerte Ritzelzähnezahl aus Niemann Getriebe allgemein, 2003, S.269, einsatzgehärtet
    i_inter = [1 2 4 8];
    z_inter1 = [32 29 25 22]; %hohe Drehzahl
    z_inter2 = [26.5 24 20.5 18]; %mittlere Drehzahl
    
    % Interpolation der Zähnezahl für Stufe 1 und 2 basierend auf Erfahrungswerten
    z(1) = round(interp1(i_inter,z_inter1,i(1),'linear'),0); %hohe Drehzahlen
    z(2) = round(interp1(i_inter,z_inter2,i(2),'linear'),0); %mittlere Drehzahlen
    
    %Für Overwrite mit bekannten Zähnezahlen aktivieren:
    % z = [23 20]; %ID.3 Zähnezahl Ritzel;
    % z = [31 24]; % Tesla Model 3 Input

else %Einstufiger Fall
    i = iges;
    if iges > 6 %Warnung wenn empfohlener Übersetzungsbereich überschritten wird
        disp("Caution, total ratio too big for one stage according to literature!")
    end
    
    % Erfahrungswerte Ritzelzähnezahl aus Niemann Getriebe allg., 2003, S.269, einsatzgehärtet
    i_inter = [1 2 4 8];
    z_inter1 = [32 29 25 22]; %hohe Drehzahl

    % Interpolation der Ritzelzähnezahl basierend auf Erfahrungswerten
    z(1) = round(interp1(i_inter,z_inter1,i(1),'linear','extrap'),0); %hohe Drehzahlen
    
    % Im einstufigen Fall nur ein Schrägungswinkel notwendig
    beta = beta(1);
end

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


%% Achsabstand a, Zahnbreite b und Normalmodul mn

M(1) = Mmax; %Drehmoment Stufe 1
M(2) = M(1)*i(1); %Drehmoment Stufe 2
for j = 1:AnzahlStufen
    %Achsabstand nach Naunheimer, Fahrzeuggetriebe, 2019, S.288 
    a(j) = 0.255 * nthroot((M(j)*1e3*(i(j)+1)^4)/i(j),3);
    
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

