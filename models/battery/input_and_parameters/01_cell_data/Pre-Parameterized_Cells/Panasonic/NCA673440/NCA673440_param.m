%% Battery parameters

%% NewParallelAssembly
NewParallelAssembly.SOC_vecCell = [0, .01, .02, .03, .04, .05, .06, .07, .08, .09, .1, .11, .12, .13, .14, .15, .16, .17, .18, .19, .2, .21, .22, .23, .24, .25, .26, .27, .28, .29, .3, .31, .32, .33, .34, .35, .36, .37, .38, .39, .4, .41, .42, .43, .44, .45, .46, .47, .48, .49, .5, .51, .52, .53, .54, .55, .56, .57, .58, .59, .6, .61, .62, .63, .64, .65, .66, .67, .68, .69, .7, .71, .72, .73, .74, .75, .76, .77, .78, .79, .8, .81, .82, .83, .84, .85, .86, .87, .88, .89, .9, .91, .92, .93, .94, .95, .96, .97, .98, .99, 1]; % Vector of state-of-charge values, SOC
NewParallelAssembly.V0_vecCell = [2.81282, 2.92952, 3.02616, 3.09882, 3.15595, 3.19736, 3.23687, 3.2696, 3.29224, 3.3165, 3.33606, 3.35376, 3.37201, 3.38673, 3.40034, 3.412, 3.42329, 3.4339, 3.44451, 3.45615, 3.466, 3.4726, 3.47919, 3.48579, 3.49239, 3.49904, 3.50604, 3.51305, 3.52005, 3.52703, 3.53325, 3.53921, 3.54518, 3.55115, 3.55712, 3.56323, 3.5696, 3.57596, 3.58233, 3.58876, 3.59531, 3.60186, 3.60839, 3.61485, 3.62133, 3.62799, 3.63464, 3.64129, 3.64862, 3.65637, 3.6643, 3.67241, 3.68051, 3.68861, 3.69647, 3.70506, 3.71368, 3.7223, 3.731, 3.74039, 3.74966, 3.75893, 3.76819, 3.77737, 3.78644, 3.79552, 3.80425, 3.81234, 3.82043, 3.82859, 3.83695, 3.84539, 3.85383, 3.86227, 3.87071, 3.87916, 3.8876, 3.89604, 3.90651, 3.91746, 3.92861, 3.94006, 3.95152, 3.96297, 3.97442, 3.98588, 3.99733, 4.0089, 4.02025, 4.03129, 4.04215, 4.0529, 4.06366, 4.07442, 4.0852, 4.09635, 4.1089, 4.1233, 4.13483, 4.15326, 4.13871]; % Open-circuit voltage, V0(SOC), V
NewParallelAssembly.V_rangeCell = [2.475, Inf]; % Terminal voltage operating range [Min Max], V
NewParallelAssembly.R0_vecCell = [.1195, .14734, .16452, .17097, .16988, .16446, .16485, .16445, .16084, .15998, .15922, .15784, .1574, .1565, .15518, .15284, .15055, .14837, .14619, .14574, .14532, .1435, .14169, .13987, .13805, .13632, .13519, .13407, .13294, .1318, .13054, .1292, .12787, .12654, .1252, .12412, .12346, .1228, .12215, .12145, .12073, .12001, .11928, .11852, .11781, .11738, .11695, .11652, .11644, .1166, .11673, .11683, .11693, .11703, .11671, .1163, .11591, .11551, .11516, .11517, .11521, .11524, .11527, .11514, .11486, .11457, .11414, .11345, .11275, .11209, .11172, .1115, .11127, .11105, .11082, .1106, .11037, .11014, .11063, .11131, .11219, .1136, .115, .11641, .11781, .11921, .12062, .12211, .12357, .1249, .12592, .12677, .12761, .12845, .12929, .13023, .13186, .13449, .1303, .12272, .05007]; % Terminal resistance, R0(SOC), Ohm
NewParallelAssembly.AHCell = 1.236; % Cell capacity, AH, A*hr

%% Additional Parameters
NewParallelAssembly.cell_type = "PrismaticGeometry";

% Source: https://industrial.panasonic.com/ww/products/pt/lithium-ion/models/NCA673440
NewParallelAssembly.cellchemistry = "NCA"; % NCA

%electrical
NewParallelAssembly.U_nom = 3.6; %nominal Voltag [V]
NewParallelAssembly.Crate_max = 4; % Assumed
%NewParallelAssembly.I_max

%physical
% Source: https://industrial.panasonic.com/ww/products/pt/lithium-ion/models/NCA673440
NewParallelAssembly.rho = 0.0203/(0.00675*0.0338*0.04035);     % Density [kg/m^3]