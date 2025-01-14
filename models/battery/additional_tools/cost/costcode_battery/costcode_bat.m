function [Cost] = costcode_bat(configs_6_BTMS_passed)

%% battery_GD inputs
battery_capacity_kWh = configs_6_BTMS_passed.SysInfo.E_sys;
h.countries = ["Germany","USA","China"];

%clear % clear memory
%clc % clear command window

%% Import Costdata
data = load("Costdata_input.mat");

%% ------------------ Extracting data from Excel sheet ------------------
% Cost of material for each cell chemistry

NMC_111 = data.NMC_111;
NMC_532 = data.NMC_532;
NMC_622 = data.NMC_622;
NMC_811 = data.NMC_811;
NCA = data.NCA;
LFP = data.LFP;
LMO = data.LMO;
Sodium = data.Sodium;

% Process Steps For Electricity cost
E_electrode_production = data.E_electrode_production;
E_cell_assembly = data.E_cell_assembly;
E_cell_finishing = data.E_cell_finishing;
E_Module_and_pack_level = data.E_Module_and_pack_level;

E_process_steps = [E_electrode_production; E_cell_assembly; E_cell_finishing; E_Module_and_pack_level];
% QUESTION = IS THERE A NEED FOR GAS COST?

% Process Steps For Labour cost
Labour_electrode_production = data.Labour_electrode_production;
Labour_cell_assembly = data.Labour_cell_assembly;
Labour_cell_finishing = data.Labour_cell_finishing;
Labour_Module_and_pack_level = data.Labour_Module_and_pack_level;

Labour_process_steps = [Labour_electrode_production; Labour_cell_assembly; Labour_cell_finishing; Labour_Module_and_pack_level];

% Country category
countries = data.countries;

Germany_Electricity_cost = data.Germany_Electricity_cost;
USA_Electricity_cost = data.USA_Electricity_cost;
China_Electricity_cost =data.China_Electricity_cost;
electricity_cost_kWh = [Germany_Electricity_cost{1, 1} USA_Electricity_cost{1, 1} China_Electricity_cost{1, 1}];
electricCostPer_Country = dictionary(countries, electricity_cost_kWh);

Germany_labour_cost = data.Germany_labour_cost;
USA_labour_cost = data.USA_labour_cost;
China_labour_cost = data.China_labour_cost;
labour_cost_kwh = [Germany_labour_cost{1,1} USA_labour_cost{1,1} China_labour_cost{1,1}];
labourCostPer_country = dictionary(countries, labour_cost_kwh);
no_of_countries = width(countries);

%Battery Capacity calculation
%battery_capacity_kWh = data.battery_capacity_kWh;

%% ---------------------- Computation Begins -------------------------------
% Arranging each cell chemistry in matrix form
%cells = {NMC_111 NMC_532 NMC_622 NMC_811 NCA LFP LMO Sodium};
cells = ["NMC_111","NMC_532","NMC_622","NMC_811","NCA","LFP","LMO","Sodium"]; %GD_main
cell_names = data.cell_names;

% Finding the no. of columns in the matrix above
cols = width(cells);

%{
fprintf('\n');
fprintf('<strong>Total Specific Battery Cost(€/kWh) Based On Country</strong>\n');
fprintf('--------------------------------------------\n');
fprintf('<strong>Battery</strong>\t\t<strong>Germany</strong>\t\t<strong>USA</strong>\t\t\t<strong>China</strong>\n');
fprintf('--------------------------------------------\n');
%}

%for c = 1:cols
for c = find(cells == configs_6_BTMS_passed.BatPara.cellchemistry) %GD_main
    % Calling the function for computation
    %fprintf('%s\t\t', cell_names{c});
    Individ_Cell_Chem_cost = cost(cells{c});

    for j=1:no_of_countries

        country_electricity_cost = electricCostPer_Country(countries(j));
        Total_Electric_cost = country_electric_cost(process_step_summation(E_process_steps, country_electricity_cost), battery_capacity_kWh);
        
        labour_cost = labourCostPer_country(countries(j));
        Total_labour_cost = country_labour_cost(process_step_summation_Labour(Labour_process_steps, labour_cost), battery_capacity_kWh);

        cost_per_country = Individ_Cell_Chem_cost + Total_Electric_cost + Total_labour_cost;
        %fprintf("%.2f\t\t", cost_per_country);
        Cost.(h.countries(j)) = cost_per_country;
    
    end
    %fprintf("\n");

end

% --------------------------------------- helpfunctions ------------------------------
%% -------------------------------------- Cell Material Cost functions ---------------
function Individ_Cell_Chem_cost = cost(cell_chemistry)
    rows = height(cell_chemistry);
    Individ_Cell_Chem_cost = 0;
    for i = 1:rows
        %fig = cell_chemistry{i,1};
        fig = cell_chemistry(i,1);
        if ~isnan(fig)
            Individ_Cell_Chem_cost = Individ_Cell_Chem_cost + fig;
        end
    end
    
end

%% -------------------------------------- Electricity Cost functions -----------------
% Function for summing each Process step
function Electric_total = process_step_summation(process_step, country_electricity_cost)
    
    %Summing each process step
    
    rows = height(process_step);
    Electric_total = 0;
    for i = 1:rows
        % Multiplying each process step by the country's electricity cost
        initial_fig = process_step{i,1}*country_electricity_cost;
        %Summation
        Electric_total = Electric_total + initial_fig;
  
    end
end

function Electric_cost = country_electric_cost(process_steps, battery_capacity_kWh) 
    %Electric_cost = process_steps*battery_capacity_kWh{1,1};
    Electric_cost = process_steps*battery_capacity_kWh; %GD_main
end

%% -------------------------------------- Labour Cost functions ----------------------
% Function for summing each Process step For Labour
function total_Labour_cost = process_step_summation_Labour(process_step, country_labour_cost)
    
    %Summing each process step
    rows = height(process_step);
    total_Labour_cost = 0;
    for i = 1:rows
        % Multiplying each process step by the country's electricity cost
        initial_fig = process_step{i,1}*country_labour_cost;
        %Summation
        total_Labour_cost = total_Labour_cost + initial_fig;
        
    end
end
%% -------------------------------------- Total Labor Cost ---------------------------
function Total_labour_cost = country_labour_cost(process_steps, battery_capacity_kWh) 
    
    %Total_labour_cost = process_steps*battery_capacity_kWh{1,1};
    Total_labour_cost = process_steps*battery_capacity_kWh; %GD_main
end

end