function T = make_sessions_table(varargin)
    P = get_parameters;
%     p=inputParser;
%     p.parse(varargin{:});
%     params=p.Results;
    T=table();
    paths = get_data_paths('session_info',true);
    for i=1:length(paths)
        load(paths{i});
        session_info.region_names = {session_info.region_names};
        T = cat(1,T,struct2table(session_info));
    end
    writetable(T,P.sessions_table_path);        
end