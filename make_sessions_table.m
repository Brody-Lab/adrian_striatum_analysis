function sessions_table = make_sessions_table(varargin)
    P = get_parameters;
%     p=inputParser;
%     p.parse(varargin{:});
%     params=p.Results;
    sessions_table=table();
    paths = get_data_paths();
    for i=1:length(paths)
        load(paths(i).session_info);
        session_info.region_names = {session_info.region_names};
        sessions_table = cat(1,sessions_table,struct2table(session_info));
    end
    writetable(sessions_table,strrep(P.sessions_table_path,'.mat','.csv'));    
    save(P.sessions_table_path,'sessions_table');        
end