function path = pathbuilder(path, prefix, suffix)
    % addapts the file pathes for GD
    arguments
        path (:,1) struct
        prefix string = ""
        suffix string = ""
    end

    names = fieldnames(path);
    for i = 1: size(names,1)
       path.(string(names(i))) = prefix + path.(string(names(i))) + suffix;
    end
end