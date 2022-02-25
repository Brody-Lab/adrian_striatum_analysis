function recordings_table = get_striatum_glm_recordings_table()
    P = get_parameters;
    recordings_table = read_recordings_log(P.recordings_path);
    recordings_table = recordings_table(recordings_table.striatum_glm==1,:);
    recordings_table.cells_file(~ismissing(recordings_table.curated_cells_file))=recordings_table.curated_cells_file(~ismissing(recordings_table.curated_cells_file));
    if any(ismissing(recordings_table.cells_file))
        warning('%g missing cells files. Skipping.\n',sum(ismissing(recordings_table.cells_file)));
    end    
    recordings_table = recordings_table(~ismissing(recordings_table.cells_file),:);
    recordings_table.cells_file = arrayfun(@(x)strrep(x,'"',''),recordings_table.cells_file);
end