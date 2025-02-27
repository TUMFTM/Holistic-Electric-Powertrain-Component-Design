%% Battery parameters

%% NewParallelAssembly
NewParallelAssembly.SOC_vecCell = [0, .01, .02, .03, .04, .05, .06, .07, .08, .09, .1, .11, .12, .13, .14, .15, .16, .17, .18, .19, .2, .21, .22, .23, .24, .25, .26, .27, .28, .29, .3, .31, .32, .33, .34, .35, .36, .37, .38, .39, .4, .41, .42, .43, .44, .45, .46, .47, .48, .49, .5, .51, .52, .53, .54, .55, .56, .57, .58, .59, .6, .61, .62, .63, .64, .65, .66, .67, .68, .69, .7, .71, .72, .73, .74, .75, .76, .77, .78, .79, .8, .81, .82, .83, .84, .85, .86, .87, .88, .89, .9, .91, .92, .93, .94, .95, .96, .97, .98, .99, 1]; % Vector of state-of-charge values, SOC
NewParallelAssembly.V0_vecCell = [2.80555, 2.85959, 2.91363, 3.00304, 3.0606, 3.09933, 3.12978, 3.18381, 3.23096, 3.26469, 3.2911, 3.31202, 3.33368, 3.35023, 3.3665, 3.38187, 3.39652, 3.40935, 3.42296, 3.43602, 3.4475, 3.45894, 3.46943, 3.47583, 3.48131, 3.48929, 3.49863, 3.50723, 3.51419, 3.52064, 3.52703, 3.53288, 3.5383, 3.54248, 3.54635, 3.55067, 3.55538, 3.56162, 3.56772, 3.57355, 3.57907, 3.58466, 3.59099, 3.59795, 3.60491, 3.61176, 3.61856, 3.62438, 3.62996, 3.63625, 3.64309, 3.64978, 3.65633, 3.66344, 3.67149, 3.67995, 3.68883, 3.69977, 3.7106, 3.72128, 3.72977, 3.73862, 3.74839, 3.75854, 3.76808, 3.77731, 3.78691, 3.79724, 3.80748, 3.817, 3.82465, 3.83227, 3.8402, 3.84937, 3.85948, 3.87001, 3.88001, 3.88909, 3.89708, 3.90554, 3.91623, 3.9276, 3.9408, 3.95458, 3.96671, 3.97597, 3.98588, 3.99921, 4.01279, 4.02654, 4.03959, 4.05214, 4.06294, 4.07257, 4.08033, 4.08863, 4.09953, 4.11328, 4.13091, 4.15201, 4.17387]; % Open-circuit voltage, V0(SOC), V
NewParallelAssembly.V_rangeCell = [2.25, Inf]; % Terminal voltage operating range [Min Max], V
NewParallelAssembly.R0_vecCell = [.17403, .14702, .12022, .10415, .08835, .07462, .06392, .06243, .06276, .06157, .06015, .0587, .05791, .05742, .05726, .05695, .05669, .0563, .05643, .05637, .05613, .05627, .0562, .0555, .05463, .05424, .05396, .05353, .05325, .053, .05273, .05239, .05185, .05104, .05016, .0496, .04917, .049, .04873, .04841, .04805, .04775, .04765, .04766, .04767, .04767, .04767, .04746, .04718, .04683, .04661, .04638, .04611, .04546, .04498, .0446, .04436, .0447, .0452, .04568, .04568, .04565, .04553, .04551, .04545, .04533, .04532, .0455, .04565, .04568, .04547, .04525, .04509, .04486, .04471, .04473, .04488, .04509, .04508, .04519, .04567, .0459, .04653, .04736, .04774, .04743, .04723, .04771, .04823, .04878, .04933, .04973, .0497, .04945, .0488, .04782, .04663, .04586, .04557, .04146, .03511]; % Terminal resistance, R0(SOC), Ohm
NewParallelAssembly.AHCell = 2.84; % Cell capacity, AH, A*hr

%% Additional Parameters
NewParallelAssembly.cell_type = "CylindricalGeometry";

% Source: https://www.akkuparts24.de/Panasonic-NCR18650PF-36V-2900mAh-Li-Ion-Zelle
NewParallelAssembly.cellchemistry = "NCA"; % NCA

%electrical
NewParallelAssembly.U_nom = 3.6; %nominal Voltag [V]
NewParallelAssembly.Crate_max = 10/2.9; % Assumed
%NewParallelAssembly.I_max

%physical
% Source: https://www.akkuparts24.de/Panasonic-NCR18650PF-36V-2900mAh-Li-Ion-Zelle
NewParallelAssembly.rho = 0.0475/((0.0186/2)^2*pi*0.0652);     % Density [kg/m^3]