function package_cells_for_tiger()
    P=get_parameters;
    recordings_table = read_recordings_log(P.recordings_path);
    curated_cells_files = recordings_table.curated_cells_file(recordings_table.striatum_glm==1);
    cells_files = recordings_table.cells_file(recordings_table.striatum_glm==1);
    cells_files(~ismissing(curated_cells_files))=curated_cells_files(~ismissing(curated_cells_files));
    fix_path = @(x)strrep(char(x),'"','');
    for i=1:length(cells_files)
        destination=fullfile(P.repository_path,'data',char(regexprep(cells_files(i),'.*Adrian(.*)','$1')));       
        if ~exist(fix_path(cells_files(i)),'file')
            error('File not found: %s.',fix_path(cells_files(i)));
        end
        if ~isdir(fileparts(destination))
            mkdir(fileparts(destination));
        end
        fprintf('Copying %s ----->\n   to %s ... \n',fix_path(cells_files(i)),fix_path(destination));tic
        copyfile(fix_path(cells_files(i)),fix_path(destination));
        fprintf('... took %s.\n-----------------\n',timestr(toc));
    end
end