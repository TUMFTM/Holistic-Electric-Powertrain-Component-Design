% costdata_evaluator v1.0
% Created for use in Global Drive 2023
% This script calculates the batterycost in €/kWh based on the data of an
% excel sheet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                              Versions                                 %
% v1.0 - creation                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           known Problems                              %
% - none                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization

clear
clc

exl = actxserver('Excel.Application');

exl.Visible = 1;    % Open Excel GUI                     

exlWkbk = exl.Workbooks;
filepath = (string(cd) + '\models\battery\additional_tools\cost\cost estimation\cost_estimation_drivetrain.xlsx');
exlFile = exlWkbk.Open(filepath);

exlSheet1= exlFile.Sheets.Item(1);

%% Cofigurate Static input
% Set the static values for the calculation here

trg.power = '200';                      %[kW]
trg.enginenr = '1';                     % As most EVs
trg.enginetype = 'PSM';       
trg.capacity = '70';                    % Set capacity to 70
trg.productionnummer = '100000';        % production of the ID3 in 2021 = 110000
trg.deprecation_period = '6';
trg.scalingfactor= 'deaktiviert';       % make sure the scalingfactor is deactivated

%% setting up field ranges
rng = struct(); % range

% possible variable values
% include new added possible variable values here
rng.enginenr_var = exlSheet1.Range('G17:G21');
rng.enginetype_var = exlSheet1.Range('H11:H12');
rng.format_var = exlSheet1.Range('H17:H19');
rng.cellchemistry_var = exlSheet1.Range('I17:I24');
rng.deprecation_period_var = exlSheet1.Range('J17:J19');
rng.country_var = exlSheet1.Range('I11:I13');
rng.scalingfactor_var = exlSheet1.Range('J11:J12');

% input ranges
rng.power = exlSheet1.Range('C12');
rng.enginenr = exlSheet1.Range('C13');
rng.enginetype = exlSheet1.Range('C14');
rng.capacity = exlSheet1.Range('C15');
rng.format = exlSheet1.Range('C16');
rng.cellchemistry = exlSheet1.Range('C17');
rng.productionnummer = exlSheet1.Range('C19');
rng.deprecation_period = exlSheet1.Range('C20');
rng.country = exlSheet1.Range('C21');
rng.scalingfactor = exlSheet1.Range('C22');

% output ranges
rng.batterycost = exlSheet1.Range('J28');

%% setting up variables
var = struct();

%static
var.enginenr_var = rng.enginenr_var.Value();
var.enginetype_var = rng.enginetype_var.Value();
var.deprecation_period_var = rng.deprecation_period_var.Value();
var.scalingfactor_var = rng.scalingfactor_var.Value();

% variated
var.cellchemistry = rng.cellchemistry_var.Value();
var.format = rng.format_var.Value();
var.country = rng.country_var.Value();

%% initialize input fields
% check out field names
fieldname.rng = fieldnames(rng);

% check if the string ends with "var"
fieldname.rng_var = {};
for i = 1:length(fieldname.rng)
    if strcmp(fieldname.rng{i}(end-2:end), 'var')
        fieldname.rng_var{end+1} = fieldname.rng{i};
    end
end

% check if value is allowend and apply it
fieldname.trg = fieldnames(trg);
for i = 1:size(fieldname.trg,1)
    if any(strcmp(string(fieldname.trg(i))+'_var',fieldname.rng_var))
        if any(strcmp(trg.(string(fieldname.trg(i))),string(var.(string(fieldname.trg(i))+'_var'))))
        rng.(string(fieldname.trg(i))).Value = trg.(string(fieldname.trg(i)));
        else
            shutdownprocedure(exlFile, exl, exlWkbk)
            error("\n!!!\nThe choosen value ''%s'' for ''%s'' is not allowed\n" + ...
                "please stay within the pool of allowed values\n!!!\n" ...
                ,string(trg.(string(fieldname.trg(i)))), string(fieldname.trg(i)));
        end
    else
        rng.(string(fieldname.trg(i))).Value = trg.(string(fieldname.trg(i)));
    end
end

%% evaluate for varied inputs

% fill struct with values
for i_format = 1:size(var.format,1)
    rng.format.Value = var.format(i_format);
    for i_cellchemistry = 1:size(var.cellchemistry,1)
        rng.cellchemistry.Value = var.cellchemistry(i_cellchemistry);
        for i_country = 1:size(var.country,1)
            rng.country.Value = var.country(i_country);

            costdata.(strfield(var.format(i_format))).(strfield(var.cellchemistry(i_cellchemistry))).(strfield(var.country(i_country)))=...
                rng.batterycost.Value()/rng.capacity.Value();   %calculates cost in [€/kWh]
        end
    end
end

% save data
save('costdata_of_celltype', 'costdata')

%% Shutdown procedure
shutdownprocedure(exlFile, exl, exlWkbk)

clearvars -except costdata

%%%%%%%%%%%%%%%%%%%%%% helpfunctions %%%%%%%%%%%%%%%%%%%%%%
%% helpfunctions that ensure allowed names of struct.fields
function struct_field_str = strfield(struct_field)
    struct_field_str = strrep(string(struct_field), '-', '_');
end

%% function that ensures a save shutdown
function shutdownprocedure(exlFile, exl, exlWkbk)
    exlFile.Close(false);       % Close workbook Excel without saving changes
    exl.Quit                    % Quit ActiveX Bridge

    % Release the ActiveX objects
    release(exlFile);
    release(exlWkbk);
    release(exl);
end