function paths=get_data_paths(varargin)
    P=get_parameters;
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('rat','',@(x)validateattributes(x,{'char','cell'},{}));
    p.addParameter('recording_name','',@(x)validateattributes(x,{'char','string'},{}));    
    p.addParameter('warn_existence',false);
    p.parse(varargin{:});
    params=p.Results;
    
    recordings_table = get_striatum_glm_recordings_table();
    idx = true(height(recordings_table),1);
    if ~isempty(params.recording_name)
        idx=recordings_table.recording_name==string(params.recording_name);
        if ~any(idx)
            error('No database match found for recording_name %s.',params.recording_name);
        end   
        if sum(idx)>1
            error('Multiple database matches found for recording_name %s.',params.recording_name);
        end          
    end
    if ~isempty(params.rat)
        idx=recordings_table.rat==string(params.rat);
        if ~any(idx)
            error('No database match found for rat %s.',params.rat);
        end
    end    
    recordings_table = recordings_table(idx,:);    
    n_sessions = height(recordings_table);
     
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