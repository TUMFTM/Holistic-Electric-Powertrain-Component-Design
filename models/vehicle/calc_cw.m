function new_cw = calc_cw(body_type, suv_flag)
%% Data origin
% Aerodynamik1 (Aerodynamik des Automobils, 5. Auflage, 
% Wolf-Heinrich Hucho, 2005)
% Site 277 (Bild 4.146)

%% Calculation
%                       A?     B       C       D           E?       F?      suv/N
    cw_table = [0.34, 0.34, 0.327, 0.314, 0.322, 0.351, 0.446];
    if suv_flag == 1
        new_cw = cw_table(7);
    else
        switch body_type
            case 'A'
                new_cw = cw_table(1);
            case 'B'
                new_cw = cw_table(2);
            case 'C'
                new_cw = cw_table(3);
            case 'D'
                new_cw = cw_table(4);
            case 'E'
                new_cw = cw_table(5);
            case 'F'
                new_cw = cw_table(6);
            case 'N'
                new_cw = cw_table(7);
        end
end