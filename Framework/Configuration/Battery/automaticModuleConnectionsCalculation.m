function configs_1_mod_all = automaticModuleConnectionsCalculation(battery_GD, input_configs)
   % Preallocate configurations array
        configs_1_mod_all = preallocate_configs_1_mod_all();

        % Step 2: Iterate through input configurations
        for ii = 1:size(input_configs, 1)
            BTMSPara = struct();
            run(input_configs{ii, 3});
            SysSpec = battery_GD.SysSpec;

            % Determine maximum serial connection size
            s_mod_raw = min(SysSpec.s_mod_max, ...
                            floor(SysSpec.U_mod_nom / battery_GD.BatPara.electrical.U_max));
            s_mod = 2 * floor(s_mod_raw / 2);

            % Determine minimum parallel connection size
            p_min_mod = max(ceil(SysSpec.I_mod_max * (1 + SysSpec.sf_parallel) / battery_GD.BatPara.electrical.I_max), ...
                            ceil(SysSpec.C_mod_min / battery_GD.BatPara.electrical.C_A));

            % Spatial arrangement of parallel connections
            epe_mod = epe_distribution(p_min_mod, SysSpec.num_higher_p_mod, ...
                                        SysSpec.mod_min_e, SysSpec.mod_max_e);

            % Iterate through configurations
            for jj = 1:size(epe_mod, 2)
                p = epe_mod(jj).p;
                s = s_mod;

                run(input_configs{ii, 2});

                SysPara = rmfield(SysPara, {'p', 's'});
                SysPara.p_mod = p;
                SysPara.s_mod = s;

                % Calculate thermal parameters
                [SysPara.thermal.transfer.A_x, SysPara.thermal.transfer.A_y, SysPara.thermal.transfer.A_z, ...
                 battery_GD.BatPara.thermal.A] = calc_cell_surface(battery_GD.BatPara);

                % Iterate through e-pe combinations
                for kk = 1:size(epe_mod(jj).pe, 2)
                    SysPara.pe_mod = epe_mod(jj).pe(kk);
                    SysPara.e_mod = epe_mod(jj).e(kk);

                    configs_1_mod_all = append_configs_step_1(battery_GD.BatPara, ...
                                                              BTMSPara, ...
                                                              SysPara, ...
                                                              SysSpec, ...
                                                              configs_1_mod_all);
                end
            end
        end
end