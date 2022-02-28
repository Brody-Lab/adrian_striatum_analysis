function cells_table = load_cells_table(varargin)
    P = get_parameters;
    load(P.cells_table_path);        
end