%% Battery parameters

%% NewParallelAssembly
NewParallelAssembly.SOC_vecCell = [0, .01, .02, .03, .04, .05, .06, .07, .08, .09, .1, .11, .12, .13, .14, .15, .16, .17, .18, .19, .2, .21, .22, .23, .24, .25, .26, .27, .28, .29, .3, .31, .32, .33, .34, .35, .36, .37, .38, .39, .4, .41, .42, .43, .44, .45, .46, .47, .48, .49, .5, .51, .52, .53, .54, .55, .56, .57, .58, .59, .6, .61, .62, .63, .64, .65, .66, .67, .68, .69, .7, .71, .72, .73, .74, .75, .76, .77, .78, .79, .8, .81, .82, .83, .84, .85, .86, .87, .88, .89, .9, .91, .92, .93, .94, .95, .96, .97, .98, .99, 1]; % Vector of state-of-charge values, SOC
NewParallelAssembly.V0_vecCell = [2.90003, 2.94341, 2.98679, 3.04987, 3.09433, 3.12791, 3.16031, 3.20369, 3.23862, 3.26689, 3.29115, 3.31383, 3.33202, 3.35112, 3.36833, 3.38413, 3.39973, 3.41137, 3.42149, 3.43136, 3.44073, 3.45027, 3.46068, 3.47113, 3.47814, 3.48494, 3.49199, 3.49905, 3.50562, 3.51187, 3.51861, 3.52535, 3.53209, 3.53818, 3.54353, 3.54849, 3.55353, 3.55855, 3.56357, 3.56831, 3.57308, 3.57787, 3.58264, 3.58755, 3.5925, 3.5985, 3.60496, 3.61149, 3.61805, 3.6246, 3.63116, 3.63855, 3.64609, 3.65336, 3.66064, 3.66791, 3.67613, 3.68522, 3.69443, 3.70367, 3.713, 3.72241, 3.73154, 3.74037, 3.7492, 3.75786, 3.76617, 3.77556, 3.78485, 3.79413, 3.8034, 3.81149, 3.81977, 3.82804, 3.83631, 3.84472, 3.85317, 3.86155, 3.86993, 3.87831, 3.88686, 3.89534, 3.90461, 3.91454, 3.92565, 3.93721, 3.94924, 3.96155, 3.97383, 3.98598, 3.99856, 4.01117, 4.02361, 4.03609, 4.04694, 4.05695, 4.06641, 4.07551, 4.08514, 4.09447, 4.10453]; % Open-circuit voltage, V0(SOC), V
NewParallelAssembly.V_rangeCell = [2.25, Inf]; % Terminal voltage operating range [Min Max], V
NewParallelAssembly.R0_vecCell = [.1741, .15823, .14313, .13285, .12044, .10641, .09605, .09162, .08758, .08349, .08021, .07738, .07487, .07338, .07189, .07117, .0704, .06898, .06755, .06631, .06509, .06388, .06326, .06281, .06182, .06095, .06025, .05955, .05877, .05807, .05761, .05714, .05668, .05611, .05546, .0549, .05439, .05388, .05337, .0528, .05223, .05165, .05106, .05046, .04986, .04945, .04913, .04883, .04854, .04826, .04797, .04782, .04765, .0473, .04695, .04661, .04638, .0463, .04625, .04621, .04624, .04628, .04628, .04623, .04617, .04599, .04573, .0457, .04568, .04565, .04563, .04535, .0451, .04486, .04461, .04446, .04435, .04423, .04411, .044, .04372, .04339, .04323, .0432, .04344, .04374, .0441, .04452, .04494, .0454, .04594, .0463, .04649, .04668, .04629, .04541, .04438, .04326, .04154, .03869, .03326]; % Terminal resistance, R0(SOC), Ohm
NewParallelAssembly.AHCell = 3.153; % Cell capacity, AH, A*hr

%% Additional Parameters
NewParallelAssembly.cell_type = "CylindricalGeometry";

% Source: https://industrial.panasonic.com/ww/products/pt/lithium-ion/models/NCR18650BD
NewParallelAssembly.cellchemistry = "NCA"; % NCA

%electrical
NewParallelAssembly.U_nom = 3.6; %nominal Voltag [V]
NewParallelAssembly.Crate_max = 4; % Assumed
%NewParallelAssembly.I_max

%physical
% Source: https://industrial.panasonic.com/ww/products/pt/lithium-ion/models/NCR18650BD
NewParallelAssembly.rho = 0.0495/((0.0185/2)^2*pi*0.0653);     % Density [kg/m^3]