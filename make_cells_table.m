function cells_table = make_cells_table(varargin)
    P = get_parameters;
%     p=inputParser;
%     p.parse(varargin{:});
%     params=p.Results;
    cells_table=table();
    paths = get_data_paths('cell_info',true);
    for i=1:length(paths)
        load(paths(i).cell_info);
        cell_info.cellno = [1:height(cell_info)]';
        cells_table = cat(1,cells_table,cell_info);
    end
    writetable(cells_table,strrep(P.cells_table_path,'.mat','.csv'));
    save(P.cells_table_path,'cells_table');        
end