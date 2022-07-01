function make_cells_agb_glm()
    
    paths = get_data_paths();
    
    for i=1:numel(paths)
        fprintf('\nFile %g.\n',i);
        Cells = load(paths(i).cells_file);
        exclude = validate_trials(Cells.Trials,'mode','agb_glm');
        Cells.Trials = remove_trials(Cells.Trials,exclude);
        Cells.spike_time_s = remove_trials_from_spikes(Cells.spike_time_s,exclude);
        new_path = strrep(paths(i).cells_file,'cells','cells_agb_glm');
        mkdir(fileparts(new_path));
        save(new_path,'-struct','Cells');        
    end



end