% -------------------------------------------------------------------------
% TU Munich - Institute for Automotive Technology (FTM)
% -------------------------------------------------------------------------
% Model for the design and analysis of a PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------

% Notes on using this script:
% (1) The input parameters can be freely defined, e.g., rated can be passed 
%     to the script, in which case the rated part below must be commented out. 
%     However, it is not mandatory to pass anything.
% (2) If the script is used without the main file, ensure the paths are set correctly.
% (3) Since no user interface for selecting the winding is invoked,
%     only the classic design can be performed (no user interaction required).

function [Entwurf, Analyse] = MEAPA_Skript(rated)

%% Select machine type
opt.Maschinentyp = rated.Maschinentyp;                                 % 'ASM', 'PMSM'

%% Input parameters for design
if(strcmp(opt.Maschinentyp,'ASM'))
    %% RATED VALUES ASM
    % % Rated power P_N [W]
    % rated.P_N   = 150000;
    % 
    % % Rated speed n_N [rpm]
    % rated.n_N   = 4300;
    % 
    % % Rated voltage U_N [V]
    % rated.U_N   = 400;
    % 
    % % Number of pole pairs p [-]
    % rated.p     = 8;
    % 
    % % Rated frequency f_N [-]
    % rated.f_N   = (rated.p * rated.n_N) / 60;
    % 
    % % Number of phases m [-]
    % rated.m     = 3;

    %% OPTIONS ASM
    % Machine type
    opt.Maschinenausfuehrung            = 'Kaefiglaeufer';                 % 'Squirrel cage rotor'
    
    % Connection
    opt.Schaltung                       = 'Star';                          % 'Star', 'Delta'
    
    % Stator winding form
    opt.Spulenform_Stator               = 'Round wire';                    % 'Round wire'
    
    % Rotor winding form
    opt.Spulenform_Rotor                = 'Round wire';                    % 'Round wire'
    
    % Stator slot shape
    opt.Nutform_Stator                  = 'Trapezoidal (angular)';         % 'Trapezoidal (angular)'
    
    % Rotor slot shape
    opt.Nutform_Rotor                   = 'Trapezoidal (angular)';         % 'Trapezoidal (angular)'
    
    % Stator cooling
    opt.Kuehlungsart                    = 'Surface cooling';               % 'Surface cooling', 'Internal or circulation cooling'
    
    % Stator iron material
    opt.Stator_Eisenmaterial.String     = 'VACOFLUX 50';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Stator_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Stator conductor material
    opt.Stator_Leitermaterial.String    = 'Copper';                        % 'Aluminum wire', 'Aluminum casting', 'Copper'
    opt.Stator_Leitermaterial = loadMaterial(opt.Stator_Leitermaterial,'Leiter');
    
    % Temperature of stator conductor material [°C]
    opt.theta_1                         = 90;
    
    % Rotor iron material
    opt.Rotor_Eisenmaterial.String      = 'VACOFLUX 50';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Rotor_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Rotor conductor material
    opt.Rotor_Leitermaterial.String     = 'Aluminum casting';              % 'Aluminum wire', 'Aluminum casting', 'Copper'
    opt.Rotor_Leitermaterial = loadMaterial(opt.Rotor_Leitermaterial,'Leiter');
    
    % Temperature of rotor conductor material [°C]
    opt.theta_2                         = 115;
    
    % Winding design mode
    opt.Mode_Wicklung                   = 'Classic';                       % 'Classic'
    
    % Winding optimization goal
    opt.Wicklungstyp                    = 'B';                             % 'A','B','C'

    %% REFERENCE VALUES ASM
    % Reference value for the relative armature length lambda [-]
    % Source: [Mueller08, p.577 - Table 9.1.3], [Meyer18, p.117]
    richt.lambda = 2.5;                                                   % between 0.6 and 1.0 for p=1, between 1.0 and 4.0 for p>1

    % Reference value for the channel width of ventilation channels l_v [m]
    % Source: [Meyer09, p.41], [Mueller08, p.585]
    richt.l_v = 0.01;                                                      % between 0.006 and 0.01

    % Reference value for the mean value of the air gap induction B_m [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_m = 0.58;                                                      % between 0.4 and 0.65

    % Reference value for current density (Stator) S_1 [A/mm^2]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    richt.S_1 = 7.0;                                                       % between 3.0 and 8.0

    % Reference value for max. allowable induction in the stator back B_1r_max [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_1r_max = 1.4;                                                  % between 1.3 and 1.65

    % Reference value for max. allowable induction in the stator teeth B_1z_max [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_1z_max = 1.8;                                                  % between 1.4 and 2.1

    % Reference value for the stator slot fill factor phi_1n [-]
    % Source: [Mueller08, p.586 - Table 9.1.6]
    richt.phi_1n = 0.5;                                                    % between 0.3 and 0.5 for round wire, between 0.35 and 0.6 for form coil or bar

    % Reference value for the winding factor (Stator) xi_1p [-]
    % Source: [Mueller08, p.596]
    richt.xi_1p = 0.96;                                                    % between 0.92 and 0.96

    % Reference value for the minimum slot pitch (Stator) tau_1n_min [m]
    % Source: [Meyer09, p.46], [Pyr14]
    richt.tau_1n_min = 0.007;                                              % between 0.007 and 0.07

    % Reference value for the iron fill factor (Stator) phi_1Fe [-]
    % Source: [Mueller08, p.599]
    richt.phi_1Fe = 0.95;                                                  % between 0.9 and 1.0

    % Reference value for the rotor current density S_2 [A/mm^2]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    richt.S_2 = 5.0;   

    % Reference value for the bar current density (Rotor) S_2s [A/mm^2]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    richt.S_2s = 4.0;                                                      % between 3.0 and 8.0 for copper, between 3.0 and 6.5 for aluminum casting, between 3.0 and 6.5 for aluminum wire

    % Reference value for the ring current density (Rotor) S_2r [A/mm^2]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    richt.S_2r = 5.0;                                                      % between 3.0 and 8.0 for copper, between 3.0 and 6.5 for aluminum casting, between 3.0 and 6.5 for aluminum wire

    % Reference value for the max. allowable induction in the rotor back B_2r_max [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_2r_max = 1.5;                                                  % between 0.4 and 1.6

    % Reference value for the max. allowable induction in the rotor teeth B_2z_max [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_2z_max = 1.9;                                                  % between 1.5 and 2.2

    % Reference value for the slot fill factor (Rotor) phi_2n [-]
    % Source: [Mueller08, p.586 - Table 9.1.6]
    % Note: Low voltage
    richt.phi_2n = 0.5;                                                    % between 0.3 and 0.5 for round wire, between 0.35 and 0.6 for form coil or bar

    % Reference value for the winding factor (Rotor) xi_2p [-]
    % Source: [Mueller08, p.596]
    % Note: No distinction between single-layer and double-layer windings,
    % as it is only required for initial estimation (calculated precisely during further design process)
    richt.xi_2p = 0.96;                                                    % between 0.92 and 0.96

    % Reference value for the minimum slot pitch (Rotor) tau_2n_min [m]
    % Source: [Meyer09, p.46], [Pyr14]
    richt.tau_2n_min = 0.007;                                              % between 0.007 and 0.07

    % Reference value for the iron fill factor (Rotor) phi_2Fe [-]
    % Source: [Mueller08, p.599]
    richt.phi_2Fe = 0.95;                                                  % between 0.9 and 1.0

elseif(strcmp(opt.Maschinentyp,'PMSM'))
    %% RATED VALUES PMSM
    % % Rated power P_N [W]
    % rated.P_N           = 115000;
    % 
    % % Rated speed n_N [rpm]
    % rated.n_N           = 4300;
    % 
    % % Rated voltage U_N [V]
    % rated.U_N           = 400;
    % 
    % % Number of pole pairs p [-]
    % rated.p             = 8;
    % 
    % % Rated frequency f_N [-]
    % rated.f_N           = (rated.p * rated.n_N) / 60;
    % 
    % % Power factor cos_phi_N [-]
    rated.cos_phi_N     = 0.9;
    % 
    % % Number of phases m [-]
    % rated.m             = 3;

    %% OPTIONS PMSM
    % Machine type
    opt.Maschinenausfuehrung            = 'SPMSM';                         % 'SPMSM', 'IPMSM (embedded)', 'IPMSM (tangential)', 'IPMSM (V-shape)';
    
    % Connection type
    opt.Schaltung                       = 'Star';                          % 'Star', 'Delta'
    
    % Stator winding form
    opt.Spulenform_Stator               = 'Round wire';                    % 'Round wire'
    
    % Stator slot shape
    opt.Nutform_Stator                  = 'Trapezoidal (angular)';         % 'Trapezoidal (angular)'
    
    % Stator cooling
    opt.Kuehlungsart                    = 'Water (direct)';                % 'Air (indirect)', 'Water (direct)'
    
    % Stator iron material
    opt.Stator_Eisenmaterial.String     = 'VACOFLUX 50';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Stator_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Stator conductor material
    opt.Stator_Leitermaterial.String    = 'Copper';                        % 'Aluminum wire', 'Aluminum casting', 'Copper'
    opt.Stator_Leitermaterial = loadMaterial(opt.Stator_Leitermaterial,'Leiter');
    
    % Temperature of stator conductor material [°C]
    opt.theta_1                         = 90;
    
    % Rotor iron material
    opt.Rotor_Eisenmaterial.String      = 'VACOFLUX 50';                   % 'M250-35A', 'M800-50A', 'VACOFLUX 48', 'VACOFLUX 50'
    opt.Rotor_Eisenmaterial = loadMaterial(opt.Stator_Eisenmaterial,'Elektroblech');
    
    % Rotor magnet material
    opt.Rotor_Magnetmaterial.String     = 'VACODYM 238 TP';                % 'VACODYM 238 TP', 'VACODYM 225 TP'
    opt.Rotor_Magnetmaterial = loadMaterial(opt.Rotor_Magnetmaterial,'Magnet');
    
    % Winding design mode
    opt.Mode_Wicklung                   = 'Classic';                       % 'Classic'
    
    % Winding optimization goal
    opt.Wicklungstyp                    = 'B';                             % 'A','B','C'

    %% REFERENCE VALUES PMSM
    % Reference value for the relative armature length lambda [-]
    % Source: [Mueller08, p.577 - Table 9.1.3], [Meyer18, p.117]
    richt.lambda = 2.5;                                                    % between 0.6 and 1.0 for p=1, between 1.0 and 4.0 for p>1

    % Reference value for the channel width of ventilation channels l_v [m]
    % Source: [Meyer09, p.41], [Mueller08, p.585]
    richt.l_v = 0.01;                                                      % between 0.006 and 0.01

    % Reference value for the amplitude of air-gap induction B_p [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_p = 0.85;                                                      % between 0.75 and 1.05

    % Reference value for the current loading A [A/mm]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    % NOTE: Either B_p or A must be specified
    %richt.A = 50.0;                                                       % between 30.0 and 120.0 for air (indirect), between 160.0 and 300.0 for water (direct)

    % Reference value for the current density (Stator) S_1 [A/mm^2]
    % Source: [Mueller08, p.580 - Table 9.1.4]
    richt.S_1 = 7.0;                                                       % between 3.0 and 7.0 for air (indirect), between 13.0 and 18.0 for water (direct)

    % Reference value for the max. allowable induction in the stator back B_1r_max [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_1r_max = 1.4;                                                  % between 1.0 and 1.5

    % Reference value for the max. allowable induction in the stator teeth B_1z_max [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_1z_max = 1.8;                                                  % between 1.6 and 2.0

    % Reference value for the stator slot fill factor phi_1n [-]
    % Source: [Mueller08, p.586 - Table 9.1.6]
    % Note: Low voltage
    richt.phi_1n = 0.5;                                                    % between 0.3 and 0.5 for round wire, between 0.35 and 0.6 for form coil or bar

    % Reference value for the winding factor (Stator) xi_1p [-]
    % Source: [Mueller08, p.596]
    % Note: No distinction between single-layer and double-layer windings,
    % as it is only required for initial estimation (calculated precisely during further design process)
    richt.xi_1p = 0.96;                                                    % between 0.92 and 0.96

    % Reference value for the minimum slot pitch (Stator) tau_1n_min [m]
    % Source: [Meyer09, p.46], [Pyr14]
    richt.tau_1n_min = 0.007;                                              % between 0.007 and 0.07

    % Reference value for the iron fill factor (Stator) phi_1Fe [-]
    % Source: [Mueller08, p.599]
    richt.phi_1Fe = 0.95;                                                  % between 0.9 and 1.0

    % Reference value for the max. allowable induction in the rotor back B_2r_max [T]
    % Source: [Mueller08, p.582 - Table 9.1.5]
    richt.B_2r_max = 1.4;                                                  % between 1.0 and 1.5

    % Reference value for the iron fill factor (Rotor) phi_2Fe [-]
    % Source: [Mueller08, p.599]
    richt.phi_2Fe = 0.95;                                                  % between 0.9 and 1.0
else
    error('Invalid input for variable "opt.Maschinentyp"');
end

%% Reassign variables
handles.rated = rated;
handles.richt = richt;
handles.opt = opt;
clear rated richt opt

%% Start design
if(strcmp(handles.opt.Maschinentyp,'ASM'))
    [handles.Entwurf] = Entwurf_ASM(handles);
elseif(strcmp(handles.opt.Maschinentyp,'PMSM'))
    [handles.Entwurf] = Entwurf_PMSM(handles);
else
    error('Invalid input for variable "Entwurf.Optionen.Maschinentyp"');
end
disp('Design completed');

%% Input parameters for analysis
    %% OPTIONS FOR LOSSES
    % Winding losses
    handles.opt.P_vw = 1;                                                  % 0, 1
    
    % Magnetic reversal losses
    handles.opt.P_vu = 1;                                                  % 0, 1
    
    % Iron loss model
    handles.opt.P_vu_Modell = 'Model approach Jordan';                     % '1'
    
    % Mechanical losses
    handles.opt.P_vme = 0;                                                 % 0, 1
    
    % Additional losses
    handles.opt.P_vzus = 1;                                                % 0, 1

    %% OPTIONS FOR CALCULATION
    % Generator operation
    handles.opt.Generator = 1;                                             % 0, 1
    
    % Max. speed [rpm]
    handles.opt.n_max = 14000;
    
    % Speed resolution
    handles.opt.n_tics = 60;
    
    % Torque resolution
    handles.opt.M_tics = 60;
    
    % Maximum voltage control [V]
    handles.opt.u_1max = handles.Entwurf.EMAG.U_1Str * sqrt(2);
    
    % Maximum current control [A]
    handles.opt.i_1max = handles.Entwurf.EMAG.I_1Str * sqrt(2);

%% Start analysis
if(strcmp(handles.Entwurf.Optionen.Maschinentyp,'ASM'))
    [handles.Analyse] = Analyse_ASM(handles);
elseif(strcmp(handles.Entwurf.Optionen.Maschinentyp,'PMSM'))
    [handles.Analyse] = Analyse_PMSM(handles);
else
    error('Invalid input for variable "handles.Entwurf.Optionen.Maschinentyp"');
end
disp('Analysis completed');

Entwurf = handles.Entwurf;
Analyse = handles.Analyse;
assignin('base','Analyse',Analyse);
end

%% Additional function
% Load material
function data = loadMaterial(var,typ)

    if(strcmp(var.String,'-'))
        data.String = var.String;
        return
    end
    switch typ
        case 'Elektroblech'
            filepath = '5_Materialien/1_Elektroblech/';
        case 'Conductor'
            filepath = '5_Materialien/2_Leiter/';
        case 'Magnet'
            filepath = '5_Materialien/3_Magnet/';
        otherwise
            error('Undefined material type')
    end

    load([filepath var.String '.mat']);
    data.String = var.String;
end