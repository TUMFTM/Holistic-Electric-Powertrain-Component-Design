% file based on transmission_ini_test.m aus Gloabl Drive

% Optimization parameters
% Name                         Range
% Transmission.gear_rat  	   6-12
% Transmission.gear_rat_S1	   2-6
% Transmission.no_mo_S1        1.5-2 mm	
% Transmission.no_mo_S2	       2-2.75 mm	
% Transmission.n_teeth_S1	   21-40
% Transmission.n_teeth_S2	   21-40

%% Inizialization of Gearbox with Gao Modell

% Design Parameter
designVector.gear_rat = 10.9;
designVector.gear_rat_S1 = 2.8;
designVector.no_mo_S1 = 1.62;
designVector.no_mo_S2 = 2.16;
designVector.n_teeth_S1 = 26;
designVector.n_teeth_S2 = 22;

% Design of Gearbox
[TransmissionDesign] = InitializationTransmission(designVector);