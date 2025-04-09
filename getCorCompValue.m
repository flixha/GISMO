function cor_comp_value = getCorCompValue(corcomp)
    % getCorCompValue returns the value of a channel / component.

    cor_comp_dict = dictionary(one="1", two="2");

    try
        cor_comp_value = cor_comp_dict(lower(corcomp));
    catch ME
        % If corcomp is not part of the cor_comp_dict, just return the input as char
        if strcmp(ME.identifier, 'MATLAB:dictionary:KeyNotFound') || ... 
                strcmp(ME.identifier, 'MATLAB:dictionary:ScalarKeyNotFound')
            cor_comp_value = char(corcomp);
        end
    end

end
