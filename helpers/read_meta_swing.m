function meta = read_meta_swing(filename)

    fid = fopen(filename,'r');
    assert(fid ~= -1, "Could not open file");
    
    meta = struct();
    
    while true
        line = fgetl(fid);
    
        if ~ischar(line) || isempty(line) || line(1) ~= '#'
            break
        end
    
        tokens = strtrim(split(line(2:end), ','));
        key = matlab.lang.makeValidName(tokens{1});
    
        values = str2double(tokens(2:end));
    
        if any(isnan(values))
            meta.(key) = tokens(2:end);
        else
            meta.(key) = values.';
        end
    end
    
    fclose(fid);

end