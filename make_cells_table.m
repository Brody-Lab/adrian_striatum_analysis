function T = make_cells_table(varargin)
    P = get_parameters;
%     p=inputParser;
%     p.parse(varargin{:});
%     params=p.Results;
    T=table();
    paths = get_data_paths('cell_info',true);
    for i=1:length(paths)
        load(paths{i});
        if i>1
            vars = intersect(T.Properties.VariableNames,cell_info.Properties.VariableNames);
            T = T(:,ismember(T.Properties.VariableNames, vars));
            cell_info = cell_info(:,ismember(cell_info.Properties.VariableNames, vars));        
        end
        T = cat(1,T,cell_info);
    end
end