function T = load_cells_table(varargin)
    P = get_parameters;
    T = readtable(P.cells_table_path);        
end