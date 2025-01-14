% Skript zur Erstellung eines Wirkungsgradsdiagramms nach MEAPA-Entwurf und
% -Analyse

% Wirkungsgradwerte aus Analyse einlesen
map = Analyse.Wirkungsgrad.eta_ges_mesh(:,2:end)*100;

%Wirkungsgradwerte kleiner 80% entfernen für bessere Darstellung
map(map<80) = NaN;

figure(5)
%Plot des Wirkungsgradkennfeldes
contour(Analyse.Betriebsdaten.n_m_mesh(60,2:end),Analyse.Betriebsdaten.M_vec,map,17)
hold on

%Plot einer x-Achse zur besseren Darstellung
plot([0 14000],[0 0],'Color',[0.5 0.5 0.5])

%Plot der Maximalleistungslinie
Mmax = max(Analyse.Betriebsdaten.M_max_mesh,[],1);
plot(Analyse.Betriebsdaten.n_m_mesh(60,:),Mmax,'k')

%Optional: Plot der Maximalleistungslinie nach Entwurfsvorgaben (rot)
% Mnenn = config_motor.power / (2*pi/60*config_motor.n_nenn);
% Mmax(Mmax>Mnenn) = Mnenn;
% plot(Analyse.Betriebsdaten.n_m_mesh(60,:),Mmax,'r')

hold off
colorbar

xlabel("Drehzahl in 1/min")
ylabel("Drehmoment in Nm")
