function T = make_cells_table(varargin)
    P = get_parameters;
%     p=inputParser;
%     p.parse(varargin{:});
%     params=p.Results;
    T=table();
    paths = get_data_paths('cell_info',true);
    for i=1:length(paths)
        load(paths{i});
        cell_info.cellno = [1:height(cell_info)]';
        T = cat(1,T,cell_info);
    end
    writetable(T,P.cells_table_path);        
end