function Trials = get_trials_from_database()
    P=get_parameters();
    
    paths = get_data_paths();
    
    parfor i=1:numel(paths)
        fprintf('Loading %s.\n',paths(i).cells_file);
        Cells = load(paths(i).cells_file);
        Trials{i} = Cells.Trials;
    end


end