function paths=get_tigress_cells_paths()
    P=get_parameters;
    recordings_table = read_recordings_log(P.recordings_path);
    curated_cells_files = recordings_table.curated_cells_file(recordings_table.striatum_glm==1);
    cells_files = recordings_table.cells_file(recordings_table.striatum_glm==1);
    cells_files(~ismissing(curated_cells_files))=curated_cells_files(~ismissing(curated_cells_files));
    fix_path = @(x)strrep(strrep(strrep(char(x),'"',''),'\',filesep),'/',filesep);
    for i=1:length(cells_files)
        paths{i}=fullfile(P.tiger_volume,'data',fix_path(char(regexprep(cells_files(i),'.*Adrian(.*)','$1'))));       
        if ~exist(path,'file')
            error('File not found: %s.',paths{i});
        end
    end
end