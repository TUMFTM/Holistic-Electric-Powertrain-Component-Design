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

secondMotor = evalin("base","secondMotor");
%% Verarbeitung für motormap2 (falls vorhanden)
if secondMotor
    % motormap2 exists and is not empty
    motormap2.skalieren = config_motor.skalieren;
    % Drehmoment- und Drehzahl-Achsen bestimmen
    nMin2 = min(min(motormap2.Speed));
    nMax2 = max(max(motormap2.Speed));
    MMin2 = min(min(motormap2.Shaft_Torque));
    MMax2 = max(max(motormap2.Shaft_Torque));
    motormap2.nVec = linspace(nMin2, nMax2, 150); % Annahme nMin2 > 0
    motormap2.MVec = [linspace(MMin2, MMax2, 120)];

    clear nMin2 nMax2 MMin2 MMax2

    % Drehmoment-Drehzahl-Gitter erstellen
    [MVecMesh2, nVecMesh2] = meshgrid(motormap2.MVec, motormap2.nVec);
    MSGID = 'MATLAB:griddata:DuplicateDataPoints';
    warning('OFF', MSGID) % Warnung für doppelte Datenpunkte deaktivieren
    % Kennfelddaten auf neue Achsen nVec und MVec fitten
    motormap2.Loss = griddata(motormap2.Speed, motormap2.Shaft_Torque, motormap2.Total_Loss, nVecMesh2, MVecMesh2, 'cubic');
    warning('ON', MSGID)

    clear nVecMesh2 MVecMesh2 MSGID

    %% Kennfeld skalieren
    if motormap2.skalieren == true
        nScale2 = config_motor.k_skal_n2; MScale2 = config_motor.k_skal_M2;
        motormap2.nVec = motormap2.nVec * nScale2; % Drehzahlachse skalieren
        motormap2.MVec = motormap2.MVec * MScale2; % Drehmomentachse skalieren
        motormap2.Loss = motormap2.Loss * nScale2 * MScale2; % Verluste skalieren
    end

end

