%% Skript zum Vorbereiten von (Motor-CAD) Kennfeldern für die Verwendung mit 
%  Simulink Lookup-Table-Blöcken + zusätzliche Skalierung n_eck und M_max

% basierend auf Matlab Beispiel "Import Efficiency Map Data from Motor-CAD"
% https://de.mathworks.com/help/sps/ug/import-efficiency-map-motorcad.html


% Informationen zur Skalierung aus Workspace laden
config_motor = evalin("base","config_motor");
motormap.skalieren = config_motor.skalieren;

% Drehmoment- und Drehzahl-Achsen bestimmen
nMin = min(min(motormap.Speed));
nMax = max(max(motormap.Speed));
MMin = min(min(motormap.Shaft_Torque));
MMax = max(max(motormap.Shaft_Torque));
motormap.nVec = linspace(nMin,nMax,150); % Annahme nMin>0
motormap.MVec = [linspace(MMin,MMax,120)];

clear nMin nMax MMin MMax


% Drehmoment-Drehzahl-Gitter erstellen
[MVecMesh,nVecMesh] = meshgrid(motormap.MVec,motormap.nVec);

MSGID = 'MATLAB:griddata:DuplicateDataPoints';
warning('OFF', MSGID) % Warnung für doppelte Datenpunkte deaktivieren
% Kennfelddaten auf neue Achsen nVec und MVec fitten
motormap.Loss = griddata(motormap.Speed,motormap.Shaft_Torque,motormap.Total_Loss,nVecMesh,MVecMesh,'cubic');
warning('ON', MSGID)

clear nVecMesh MVecMesh MSGID

%% Kennfeld skalieren
if motormap.skalieren == true
    nScale = config_motor.k_skal_n; MScale = config_motor.k_skal_M;
    motormap.nVec = motormap.nVec*nScale; %Drehzahlachse skalieren
    motormap.tVec = motormap.MVec*MScale; %Drehmomentachse skalieren
    motormap.Loss = motormap.Loss*nScale*MScale; %Verluste skalieren
end