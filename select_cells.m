function cells_table = select_cells(varargin)
    
    %% parse and validate inputs
    p=inputParser;
    p.addParameter('DV',[-Inf Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2}));
    p.addParameter('AP',[-Inf Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2}));
    p.addParameter('ML',[0 Inf],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2,'>=',0}));
    p.addParameter('regions',{},@(x)validateattributes(x,{'char','cell'},{}));
    p.addParameter('spike_width_ms',[],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2,'>=',0}));
    p.addParameter('reliability',[],@(x)validateattributes(x,{'numeric'},{'increasing','numel',2}));
    p.addParameter('D2Phototagging',[],@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('is_in_dorsal_striatum',[],@(x)validateattributes(x,{'logical'},{'scalar'}));        
    p.addParameter('recording_name',"",@(x)validateattributes(x,{'string'},{'scalar'}));
    p.parse(varargin{:});
    select_cells_params=p.Results;  

    %% make a cells_table with sessions columns
    cells_table = load_cells_table();
    sessions_table = load_sessions_table();
    cells_table = join(cells_table,sessions_table);
    
    %% 
    cells_table = validate_ranges(cells_table,select_cells_params);

end


function cells_table = validate_ranges(cells_table,params)
    
    select_cells_fields = fieldnames(params);
    cells_table_fields = cells_table.Properties.VariableNames;
    include = true(height(cells_table),length(select_cells_fields));                
    for f=1:length(select_cells_fields)       
        prop = select_cells_fields{f};
        val = params.(select_cells_fields{f});        
        if ismember(prop,cells_table_fields)
            if isempty(val)
            elseif isstring(val)
                include(:,f) = cells_table.(prop) == val;
            elseif isscalar(val) && ~iscell(val)
                
                if isnan(val)
                    include(:,f) = isnan(cells_table.(prop));
                elseif islogical(val)
                    if ~islogical(cells_table.(prop))
                        warning('select_cells: property %s of cells_table is not logical but selection criteria is. Check results carefully.',prop);
                    end
                    include(:,f) = cells_table.(prop) == double(val) ;
                else
                    include(:,f) = cells_table.(prop) == val;
                end
                
            elseif isnumeric(val) && numel(val)==2
                
                include(:,f) = cells_table.(prop)>val(1) & cells_table.(prop)<=val(2); % bounds are [ )
            elseif iscell(val)
                
                error(' haven''t implemented this yet.');
            else 
                warning('select_cells: Value for %s has an unrecognized type or size.',val);
            end
        else
            error('%s is not a column of either the cells_table or sessions_table.',prop);
        end
        if any(~include(:,f))
            fprintf('select_cells:Removed %g cells due to %s criterion.\n',sum(~include(:,f)),prop);
        end
    end   
    include = all(include,2);
    if ~any(include)
        warning('select_cells: No cells matched criteria.');
    elseif all(include)
        fprintf('select_cells: ALL %g cells in database matched criteria.\n',sum(include));
    else
        fprintf('select_cells: %g of %g cells in database (%g%%) matched criteria.\n',sum(include),numel(include),round(100*mean(include)));
    end
    cells_table = cells_table(include,:);
end
