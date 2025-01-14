% Script for gear data calculations
% Sources: Niemann, Getriebe allgemein, 2003, S.58 und 276f.
%          FZG, Formelsammlung zum Modul Maschinenelemente, 2020, S.102ff.
%          Gruss, Berechnungsgrundlagen Stirnradverzahnungen (online verfügbar)

alphawt = atand(tand(alpha)/cosd(beta)); % Degrees: Operating pressure angle
betab = atand(tand(beta)*cosd(alphawt)); % Degrees: Helix angle at the base circle
mt = mn/cosd(beta); % mm: Transverse module

u = z2/z1; % --: Gear ratio; transmission ratio

d1 = mt*z1; % mm: Pitch diameter of gear 1
d2 = mt*z2; % mm: Pitch diameter of gear 2

dw1 = d1; % mm: Rolling diameter of gear 1
dw2 = d2; % mm: Rolling diameter of gear 2

db1 = d1 * cosd(alphawt); % mm: Base circle diameter of gear 1
db2 = d2 * cosd(alphawt); % mm: Base circle diameter of gear 2

a = (d1 + d2)/2; % mm: Center distance

haP = mn; % mm: Addendum height
da1 = d1 + 2*haP; % mm: Tip diameter of gear 1
da2 = d2 + 2*haP; % mm: Tip diameter of gear 2

galpha = 0.5 * (sqrt(da1^2 - db1^2) + sqrt(da2^2 - db2^2)) - a * sind(alphawt); % mm: Path of contact
galphaa1 = 0.5 * db1 * (((da1/db1)^2 - 1)^0.5 - tand(alphawt)); % mm: Contact path on pinion tip
eps_alpha = galpha / (mt * pi * cosd(alphawt)); % Total contact ratio
eps_1 = galphaa1 / (mt * pi * cosd(alphawt)); % Partial contact ratio (pinion tip)
eps_2 = eps_alpha - eps_1; % Partial contact ratio (gear tip)
