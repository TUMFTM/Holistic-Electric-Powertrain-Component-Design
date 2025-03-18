function [m_gearbox, V_gearbox] = calcgearboxweight(GB_Topology, GB, config_gearbox)

 if strcmp(GB_Topology, 'achsparallel')
        % from: Application of the TOPSIS Method for MultiObjective Optimization of a Two-Stage Helical Gearbox
        % the volume calculations are reused for the koaxial case for
        % estimation
       
        % Gearbox housing dimensions
        d1_stage1 = GB.d1(1); % Diameter of pinion (Stage 1)
        d2_stage1 = GB.d2(1); % Diameter of wheel (Stage 1)
        b1_stage1 = GB.b(1); % Tooth width (Stage 1)
        d1_stage2 = GB.d1(2); % Diameter of pinion (Stage 2)
        d2_stage2 = GB.d2(2); % Diameter of wheel (Stage 2)
        b1_stage2 = GB.b(2); % Tooth width (Stage 2)
        d_shaft_1 = GB.d(1); % Shaft diameter (Stage 1)
        d_shaft_2 = GB.d(2); % Shaft diameter (Stage 2)
        d_shaft_3 = GB.d(3); % Shaft diameter (Stage 3)
    
        % Constants
        rho_g = 7850; % Density of gears (kg/m^3)
        rho_steel = 7850; % Density of steel (kg/m^3)
        rho_alualloy = 2700; % Density of aluminum alloy (kg/m^3)
    
        % Gearbox housing dimensions
        L = (d1_stage1 + d2_stage1 / 2 + d1_stage2 / 2 + d2_stage2 / 2 + 22.5) / 0.975;
        S_G = 0.005 * L + 4.5;
        H = max(d2_stage1, d2_stage2) + 8.5 * S_G;
        B = b1_stage1 + b1_stage2 + 6 * S_G;
         
        % Shaft lengths
        L_1 = B + 1.2 * d_shaft_1;
        L_2 = B;
        L_3 = B + 1.2 * d_shaft_3;
    
        % Volumes of shafts
        V_shaft_1 = (pi / 4) * d_shaft_1^2 * L_1;
        V_shaft_2 = (pi / 4) * d_shaft_2^2 * L_2;
        V_shaft_3 = (pi / 4) * d_shaft_3^2 * L_3;
        V_shafts = (V_shaft_1 + V_shaft_2 + V_shaft_3)/1000000000;
    
        % Volumes of gears
        V_gear_1_stage1 = (pi / 4) * (d1_stage1^2 - d_shaft_1^2) * b1_stage1;
        V_gear_2_stage1 = (pi / 4) * (d2_stage1^2 - d_shaft_2^2) * b1_stage1;
        V_gear_1_stage2 = (pi / 4) * (d1_stage2^2 - d_shaft_2^2) * b1_stage2;
        V_gear_2_stage2 = (pi / 4) * (d2_stage2^2 - d_shaft_3^2) * b1_stage2;
        V_gears = (V_gear_1_stage1 + V_gear_2_stage1 + V_gear_1_stage2 + V_gear_2_stage2)/1000000000;
    
        % Housing volumes
        V_gearbox = L * H * B/1000000000;
        V_A = L * H * S_G/1000000000;
        V_B = L * b1_stage1 * 1.5 * S_G/1000000000;
        V_C = (B - 2 * S_G) * H * S_G/1000000000;
    
        % Mass calculations
        e_1 = 1; e_2 = 0.6;
        m_g1 = rho_g * ((pi * e_1 * d1_stage1^2 * b1_stage1) / 4 + (pi * e_2 * d2_stage1^2 * b1_stage1) / 4)/1000000000;
        m_g2 = rho_g * ((pi * e_1 * d1_stage2^2 * b1_stage2) / 4 + (pi * e_2 * d2_stage2^2 * b1_stage2) / 4)/1000000000;
        m_shafts = rho_steel * V_shafts;
        m_housing = rho_alualloy * 2 * (V_A + V_B + V_C);
        m_gears = m_g1 + m_g2;
        m_gearbox = m_shafts + m_housing + m_gears;
        
     
    elseif strcmp(GB_Topology, 'koaxial')
        m_gearbox = 0.199*((GB.iges*config_gearbox.M_max)^0.669*GB.Stufen^0.334);    % Source: Fahrzeuggetriebe, Naunheimer 3 Auflage Seite 55
        

        % Gearbox housing dimensions
        d1_stage1 = GB.d1(1); % Diameter of pinion (Stage 1)
        d2_stage1 = GB.d2(1); % Diameter of wheel (Stage 1)
        b1_stage1 = GB.b(1); % Tooth width (Stage 1)
        d1_stage2 = GB.d1(2); % Diameter of pinion (Stage 2)
        d2_stage2 = GB.d2(2); % Diameter of wheel (Stage 2)
        b1_stage2 = GB.b(2); % Tooth width (Stage 2)
        d_shaft_1 = GB.d(1); % Shaft diameter (Stage 1)
        d_shaft_2 = GB.d(2); % Shaft diameter (Stage 2)
        d_shaft_3 = GB.d(3); % Shaft diameter (Stage 3)
    
        % Constants
        rho_g = 7850; % Density of gears (kg/m^3)
        rho_steel = 7850; % Density of steel (kg/m^3)
        rho_alualloy = 2700; % Density of aluminum alloy (kg/m^3)
    
        % Gearbox housing dimensions
        L = (d1_stage1 + d2_stage1 / 2 + d1_stage2 / 2 + d2_stage2 / 2 + 22.5) / 0.975;
        S_G = 0.005 * L + 4.5;
        H = max(d2_stage1, d2_stage2) + 8.5 * S_G;
        B = b1_stage1 + b1_stage2 + 6 * S_G;

        V_gearbox = L * H * B/1000000000;
        
     
    else
        m_gearbox = 10; %Placeholder for equation
        V_gearbox = 0; 
   
end
end