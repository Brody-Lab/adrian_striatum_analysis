function T = load_sessions_table(varargin)
    P = get_parameters;
    T = readtable(P.sessions_table_path);        
end