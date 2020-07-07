function [saved_cells,file_paths] = get_saved_cellnos(path)
    % returns the saved cellnos for a run plus the paths to the stats
    % structures
    % path should be the path to the fit run (i.e. the parent folder of "stats"
    stats_dir = fullfile(path,'stats');
    if isdir(stats_dir)
        tmp=dir(stats_dir);
    else
        error('No "stats" folder in that path.');
    end
    saved_cells = cellfun(@(x)str2num(x),regexprep({tmp.name},'.*cell(.*).mat','$1'),'uni',0);
    is_a_cell_file = ~cellfun(@(x)isempty(x),saved_cells);
    tmp=tmp(is_a_cell_file);
    saved_cells = [saved_cells{is_a_cell_file}];
    file_paths={tmp.name};
    file_paths = cellfun(@(x)fullfile(stats_dir,x),file_paths,'uni',0);
    [saved_cells,sort_idx] = sort(saved_cells);
    file_paths=file_paths(sort_idx);
end
