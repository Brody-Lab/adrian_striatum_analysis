function paths=get_data_paths(varargin)
    P=get_parameters;
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('rat','',@(x)validateattributes(x,{'char','cell'},{}));
    p.addParameter('recording_name','',@(x)validateattributes(x,{'char','string'},{}));    
    p.addParameter('warn_existence',false);
    p.addParameter('joined',false);    
    p.addParameter('data_path',P.data_path,@(x)validateattributes(x,{'char','string'},{'nonempty'}));
    p.parse(varargin{:});
    params=p.Results;
    P.data_path = params.data_path;
    if ~isfolder(P.data_path)
        error('Main data path %s does not exist.',P.data_path);
    end
    recordings_table = get_striatum_glm_recordings_table();
    idx = true(height(recordings_table),1);
    joined=false;
    if ~isempty(params.recording_name)
        [rat,date,suffix] = parse_recording_name(recordings_table.recording_name);        
        [this_rat,this_date,this_suffix] = parse_recording_name(params.recording_name);
        if isempty(this_suffix)
            idx=rat==this_rat&date==this_date;
        else
            idx=rat==this_rat&date==this_date&suffix==this_suffix;
        end            
        if ~any(idx)
            error('No database match found for recording_name %s.',params.recording_name);
        end   
        if sum(idx)>1
            joined=true;
        end  
        if ~isempty(params.rat)
            error('Specify only one of "rat" or "recording_name".');
        end
    elseif ~isempty(params.rat)
        idx=recordings_table.rat==string(params.rat);
        if ~any(idx)
            error('No database match found for rat %s.',params.rat);
        end
    end    
    recordings_table = recordings_table(idx,:);    
    n_sessions = height(recordings_table);
    if joined
        paths.cells_file = fullfile(P.data_path,'cells','joined',params.recording_name,[char(params.recording_name),'_Cells.mat']);
        for i=n_sessions:-1:1
            paths.original_cells_file{i} = recordings_table.cells_file(i);
        end
        paths.recording_name = params.recording_name;
        paths.date = recordings_table.date(i);  
        paths.parent_dir = fileparts(paths.cells_file); 
        paths.fit_path = fullfile(P.fit_path,'fits','joined',params.recording_name);
        paths.all_exist=true;
        if ~exist(paths.cells_file,'file')
            paths.all_exist = false;
            if params.warn_existence
                warning('Expected file in local database does not exist: %s.',paths.cells_file);
            end
        end           
    else
        for i=n_sessions:-1:1
            paths(i).cells_file=fullfile(P.data_path,'cells',recordings_table.recording_name(i),[char(recordings_table.recording_name(i)),'_Cells.mat']);    
            if ~exist(paths(i).cells_file,'file')
                error('');
            end
            paths(i).original_cells_file = recordings_table.cells_file(i);
            paths(i).rat_name = recordings_table.rat_name(i);
            paths(i).recording_name = recordings_table.recording_name(i);
            paths(i).date = recordings_table.date(i);  
            paths(i).parent_dir = fileparts(paths(i).cells_file);
            paths(i).cell_info = fullfile(paths(i).parent_dir,'cell_info.mat');
            paths(i).session_info = fullfile(paths(i).parent_dir,'session_info.mat');   
            paths(i).fit_path = fullfile(P.fit_path,'fits',recordings_table.recording_name(i));
            paths(i).all_exist=true;
            if ~exist(paths(i).cells_file,'file')
                paths(i).all_exist = false;
                if params.warn_existence
                    warning('Expected file in local database does not exist: %s.',paths(i).cells_file);
                end
            end   
            if ~exist(paths(i).cell_info,'file')
                paths(i).all_exist = false;        
                if params.warn_existence            
                    warning('Expected file in local database does not exist: %s.',paths(i).cell_info);
                end
            end   
            if ~exist(paths(i).session_info,'file')
                paths(i).all_exist = false;       
                if params.warn_existence
                    warning('Expected file in local database does not exist: %s.',paths(i).session_info);
                end
            end           
        end
    end
end