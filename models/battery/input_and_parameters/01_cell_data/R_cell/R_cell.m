%% Disclaimer

% This is no real cell data! Use at your own risk!

% The assumptions for cell capacities and dimensions are based on the
% following publication: Fraunhofer-Allianz Batterien (Hg.) (2017): 
% Enwicklungsperspektiven für Zellformate von Lithium-Ionen-Batterien 
% in der Elektromobilität.

%% Cell type

BatPara.name = mfilename;       % Get cell name from filename (for better overview)

BatPara.cell_type = 'Cyl';      % Select either 'Pouch', 'Pris' or 'Cyl'

BatPara.cellchemistry = "NMC_811"; % [NMC_111, NMC_532, NMC_622, NMC_811, NCA, LFP, LMO, Sodium ]

%% Electrical parameters (edit here)

% Specify Electrical LIB-Parameters. 

% Specify the BatPara structure. Details can be found at https://github.com/TUMFTM/sim_battery_system

% For reasons of simplicity and comparability we assume the dynamical values 
% used in the initial release of the electrical system simulation which can 
% be found at https://github.com/TUMFTM/sim_battery_system.

run cell_parameters_NCR18650PF


t_cell_capacity = 4.900;   % Cell capacigy C_A  in Ah

% Scale sample parameter set to approximate dynamic behavior
[BatPara] = scale_BatPara(BatPara, t_cell_capacity);

BatPara.electrical.C_A = 4.900;                     % Cell capacigy C_A  in Ah
BatPara.electrical.U_max = 4.2;                     % Maximum allowed cell voltage in V IN THE LATER USE-CASE
BatPara.electrical.U_min = 2.5;                     % Minimum allowed cell voltage in V IN THE LATER USE-CASE
BatPara.electrical.U_nom = 3.63;                    % Nominal cell voltage in V
BatPara.electrical.Crate_max = 3;                   % maximal allowed C-Rate "agressive cycling"
BatPara.electrical.I_max = BatPara.electrical.Crate_max * t_cell_capacity;     % Maximum cell current in A (assumption: 2C) --> Warning, this is very high for current gen cells!
%Source: https://www.ancoo-battery.com/en/product/INR21700-50G.html

%% Physical and thermal parameters (edit here)

% Note: Because we use a different thermal simulation model in this model
% compared to https://github.com/TUMFTM/sim_battery_system, we don't use
% the thermal parameters and the format of this model, but our own.

% Source: Xia, Et al. (2018): A reliability design method for a lithium-ion battery pack considering the thermal disequilibrium in electric vehicles. 
% In: Journal of Power Sources 386, S. 10–20. DOI: 10.1016/j.jpowsour.2018.03.036.

BatPara.physical.c = 696.07;        % specific heat capacity [J/(kg*K)]

% Cell dimensions. Note: For cylindrical cells two of those values must be
% the same --> This will be considered as diameter. The diameter must be
% identical or smaller as the length!

BatPara.physical.dim_x = 21.1e-3;     % Cell dimension in x-direction
BatPara.physical.dim_y = 70.15e-3;    % Cell dimension in y-direction
BatPara.physical.dim_z = 21.1e-3;     % Cell dimension in z-direction

% Cell mass (Calculated from the cell's dimensions and density)
BatPara.physical.m = 69.5e-3;
BatPara.physical.rho = get_cell_density(BatPara.physical.dim_x, BatPara.physical.dim_y, BatPara.physical.dim_z, BatPara.physical.m, BatPara.cell_type);     % Density [kg/m^3]

% new Simscape Parameters

%% Clear temporary vars

clearvars t_*