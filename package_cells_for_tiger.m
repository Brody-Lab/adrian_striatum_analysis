function package_cells_for_tiger(varargin)
    P=get_parameters;
    p=inputParser;
    p.addParameter('make_cell_info',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('remake',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.parse(varargin{:});
    params = p.Results;
    recordings_table = read_recordings_log(P.recordings_path);
    tagged = recordings_table.D2Phototagging==1;
    tagged = tagged(recordings_table.striatum_glm==1);
    laser_power_mW = recordings_table.laser_power_mW;
    laser_power_mW = laser_power_mW(recordings_table.striatum_glm==1);    
    curated_cells_files = recordings_table.curated_cells_file(recordings_table.striatum_glm==1);
    cells_files = recordings_table.cells_file(recordings_table.striatum_glm==1);
    cells_files(~ismissing(curated_cells_files))=curated_cells_files(~ismissing(curated_cells_files));
    if any(ismissing(cells_files))
        warning('%g missing cells files. Skipping.\n',sum(ismissing(cells_files)));
    end    
    tagged = tagged(~ismissing(cells_files));
    laser_power_mW = laser_power_mW(~ismissing(cells_files));
    cells_files = cells_files(~ismissing(cells_files));
    fix_path = @(x)strrep(char(x),'"','');
    for i=1:length(cells_files)
%         if ~tagged(i)
%             continue
%         end
        destination=fullfile(P.data_path,'cells',char(regexprep(cells_files(i),'.*Adrian(.*)','$1')));       
        if ~exist(fix_path(cells_files(i)),'file')
            error('File not found: %s.',fix_path(cells_files(i)));
        end
        if ~isdir(fileparts(destination))
            mkdir(fileparts(destination));
        end
        if isfile(fix_path(destination)) && ~params.remake
            fprintf('%s exists.\n',fix_path(destination));         
        else
            fprintf('Copying %s ----->\n   to %s ...',fix_path(cells_files(i)),fix_path(destination));tic
            copyfile(fix_path(cells_files(i)),fix_path(destination));
            fprintf(' took %s.\n',timestr(toc));        
        end
        if params.make_cell_info
            [parent,~,~] = fileparts(fix_path(destination));
            cell_info_path = fullfile(parent,'cell_info.mat');            
            if isfile(cell_info_path) && ~params.remake 
                fprintf('%s exists.\n-----------------\n',cell_info_path);                         
            else
                fprintf('Loading Cells file and making cell_info structure ...');tic;
                Cells = load(fix_path(destination));
                fields=fieldnames(Cells);
                if length(fields)==1
                    Cells=Cells.(fields{1}); % if not saved with -struct flag
                end   
                Cells.laser_power_mW = laser_power_mW(i);
                cell_info = make_cell_info(Cells,tagged(i));
                save(cell_info_path,'cell_info');
                fprintf(' took %s.\n-----------------\n',timestr(toc));
            end
        end
    end
end