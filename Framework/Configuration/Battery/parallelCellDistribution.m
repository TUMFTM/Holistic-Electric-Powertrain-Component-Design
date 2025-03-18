function configs_1_mod_all = parallelCellDistribution(battery_GD, input_configs)
    
    % parallelCellDistribution
    % Generates all possible configurations for the parallel arrangement of cells
    % in a battery module, taking into account spatial and thermal parameters.
    %
    % Inputs:
    %   battery_GD: Struct containing battery parameters, system specifications,
    %               and system-specific parameters required for calculations.
    %
    % Outputs:
    %   configs_1_mod_all: Cell array containing all possible configurations for
    %                      the parallel arrangement of cells in the module.

    % Extract battery parameters
    BatPara = battery_GD.BatPara;

    % Preallocate the cell array for configurations
    configs_1_mod_all = preallocate_configs_1_mod_all();

    % Step 1: Spatial arrangement of the parallel connection
    % The cells can be arranged in a 1D grid (s*p) or across all three dimensions.
    % The 'e*pe' connection implies that p = e * pe, where:
    %   - e: Number of cells in one dimension
    %   - pe: Number of cells in the other dimension
    % Total number of cells in a module: n = e * pe * s.
    % Check if the field 'p_mod' exists in SysPara

    % Use p_mod_value in your epe_distribution function
    epe_mod = epe_distribution(battery_GD.SysPara.p_mod, ...
                               battery_GD.SysSpec.num_higher_p_mod, ...
                               battery_GD.SysSpec.mod_min_e, ...
                               battery_GD.SysSpec.mod_max_e);
   
  
    run(input_configs{1,3});
     
    % Step 2: Iterate through all possible configurations of e and pe
    for jj = 1:size(epe_mod, 2)
        
        p=battery_GD.SysPara.p_mod;
        s=battery_GD.SysPara.s_mod;
        run(input_configs{1,2});

        SysPara = battery_GD.SysPara;
        
        % Step 2.1: Calculate additional thermal cell parameters
        [SysPara.thermal.transfer.A_x, SysPara.thermal.transfer.A_y, SysPara.thermal.transfer.A_z, ...
         BatPara.thermal.A] = calc_cell_surface(BatPara);

        % Step 2.2: Iterate through all e-pe combinations for the current p
        for kk = 1:size(epe_mod(jj).pe, 2)
            
            % Update system parameters for the current configuration
            SysPara.pe_mod = epe_mod(jj).pe(kk);
            SysPara.e_mod = epe_mod(jj).e(kk);
              
            % Save the configuration to the cell array
            configs_1_mod_all = append_configs_step_1(BatPara, ...
                                                      BTMSPara, ...
                                                      SysPara, ...
                                                      battery_GD.SysSpec, ...
                                                      configs_1_mod_all);
        end
    end
end