function path = correct_file_path(path)
    P=get_parameters;
    if P.on_tiger
        path  = strrep(path,P.pc_data_path,P.tiger_data_path);
    else
        path  = strrep(path,P.tiger_data_path,P.pc_data_path);        
    end
    path = fixfilesep(path);
end