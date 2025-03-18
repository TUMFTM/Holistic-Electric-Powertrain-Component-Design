function config_x = read_json_config(file_path)
    % Read json formatted table of config data and import it to Matlab workspace
    fid = fopen(file_path, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    config_x = jsondecode(str);
end