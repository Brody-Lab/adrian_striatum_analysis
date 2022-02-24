function paths=get_data_paths(varargin)
    P=get_parameters;
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('data_path',fullfile(P.data_path,'cells'));
    p.addParameter('parent_dir',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('rat','',@(x)validateattributes(x,{'char','cell'},{}));
    p.addParameter('must_exist',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('cell_info',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('session_info',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.parse(varargin{:});
    params=p.Results;
    recordings_table = read_recordings_log(P.recordings_path);
    curated_cells_files = recordings_table.curated_cells_file(recordings_table.striatum_glm==1);
    cells_files = recordings_table.cells_file(recordings_table.striatum_glm==1);
    rat_names = recordings_table.rat_name(recordings_table.striatum_glm==1);
    cells_files(~ismissing(curated_cells_files))=curated_cells_files(~ismissing(curated_cells_files));
    if any(ismissing(cells_files))
        warning('%g missing cells files. Skipping.\n',sum(ismissing(cells_files)));
    end
    cells_files = cells_files(~ismissing(cells_files));
    rat_names = rat_names(~ismissing(cells_files));
    fix_path = @(x)strrep(strrep(strrep(char(x),'"',''),'\',filesep),'/',filesep);
    for i=1:length(cells_files)
        paths{i}=fullfile(params.data_path,fix_path(char(regexprep(cells_files(i),'.*Adrian(.*)','$1'))));   
        parent = fileparts(paths{i});
        if params.parent_dir
            paths{i} = parent;
            if ~isdir(paths{i}) && params.must_exist
                error('Directory not found: %s.',paths{i});                  
            end            
        else
            if params.cell_info
                paths{i} = fullfile(parent,'cell_info.mat');
            elseif params.session_info
                paths{i} = fullfile(parent,'session_info.mat');
            end
            if ~exist(paths{i},'file') && params.must_exist
                error('File not found: %s.',paths{i});
            end
        end
    end
    if ~isempty(params.rat)
        if ischar(params.rat)
            params.rat={params.rat};
        end
        paths = paths(ismember(rat_names,params.rat));
    end
end