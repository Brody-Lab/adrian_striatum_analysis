function recordings_table = get_striatum_glm_recordings_table()
    P = get_parameters;
    recordings_table = read_recordings_log(P.recordings_path);
    recordings_table = recordings_table(recordings_table.striatum_glm==1,:);
    use_curated_if_available_idx = ~ismissing(recordings_table.curated_cells_file) & ~recordings_table.UberPhys;    % I curated some of the uberphys sessions but they had a bug in them and were later resorted
    recordings_table.cells_file(use_curated_if_available_idx)=recordings_table.curated_cells_file(use_curated_if_available_idx);
    if any(ismissing(recordings_table.cells_file))
        warning('%g missing cells files. Skipping.\n',sum(ismissing(recordings_table.cells_file)));
    end    
    recordings_table = recordings_table(~ismissing(recordings_table.cells_file),:);
    recordings_table.cells_file = arrayfun(@(x)strrep(x,'"',''),recordings_table.cells_file);
    recordings_table.cells_file = strtrim(recordings_table.cells_file);
end