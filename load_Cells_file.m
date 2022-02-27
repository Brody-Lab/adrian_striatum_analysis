function Cells = load_Cells_file(recording_name,use_local)
    if nargin==1
        use_local = true;
    end       
    paths = get_data_paths('recording_name',recording_name);
    if use_local
        if isfile(paths.cells_file)
            fprintf('Loading Cells file from local copy @ %s ... \n',paths.cells_file);tic;
            Cells = load(paths.cells_file);
            fprintf(' took %s.\n-----------------\n',timestr(toc));      
        else
            fprintf('Local file does not exist yet. Loading Cells file from original location @ %s... \n',paths.original_cells_file);tic;
            Cells = load(paths.original_cells_file);
            fprintf('   ... took %s.\n-----------------\n',timestr(toc));   
        end
    else
        fprintf('Loading Cells file from original location @ %s... \n',paths.original_cells_file);tic;
        Cells = load(paths.original_cells_file);
        fprintf('   ... took %s.\n-----------------\n',timestr(toc));                 
    end 
    if isfield(Cells,'Cells')
        Cells=Cells.Cells;
    end
end