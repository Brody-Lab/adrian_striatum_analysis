function sessions_table = load_sessions_table(varargin)
    P = get_parameters;
    load(P.sessions_table_path);        
end