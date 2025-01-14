function struct_1 = config_enforcer(struct_1, struct_2)
%function for ensuring the usage of the the GD config
    overlap = intersect(fieldnames(struct_1), fieldnames(struct_2)); 
    for i = 1:size(overlap,1)
        struct_1.(string(overlap(i)))=struct_2.(string(overlap(i)));
    end
    
    struct_1.dim_x_mod_max = struct_2.dim_x_sys_max*(10/11);
    struct_1.dim_y_mod_max = struct_2.dim_y_sys_max*(10/11);
    struct_1.dim_z_mod_max = struct_2.dim_z_sys_max*(10/11);
end