function include = apply_inclusion_criteria(S,varargin)
    
%apply_inclusion_criteria Determine which cells in a Cells file satisfy the striatum analysis inclusion criteria
%   INCLUDE = APPLY_INCLUSION_CRITERIA(cells_table) returns a logical array 
%   with values of 1 for units in the cells table satisfying all inclusion criteria, 
%   and values of 0 otherwise.
%
%   INCLUDE = APPLY_INCLUSION_CRITERIA(CELLS,'MIN_VALUE_FIELD',MIN_VALUE_FIELD_NAME)
%   uses the column MIN_VALUE_FIELD_NAME in the
%   unit_inclusion_critera.csv table for the mininum allowed
%   parameter values instead of the default column.
%
%   INCLUDE = APPLY_INCLUSION_CRITERIA(CELLS,'MAX_VALUE_FIELD',MAX_VALUE_FIELD_NAME)
%   uses the column MAX_VALUE_FIELD_NAME in the
%   unit_inclusion_critera.csv table for the maxinum allowed
%   parameter values instead of the default column.
%
%   Examples:
%
%   >> S = load_cells_table();
%   >> apply_inclusion_criteria(S)
%   ans = 
% 
%       18762×1 logical array
% 
%           0
%           1
%           1
%           1
%           
%           ...

    %% function argument parsing and validation
    T = get_inclusion_criteria();
    p=inputParser;
    p.addParameter('min_value_field','min_value',@(x)validateattributes(x,{'string','char'},{'nonempty'}));
    p.addParameter('max_value_field','max_value',@(x)validateattributes(x,{'string','char'},{'nonempty'}));  
    p.parse(varargin{:});
    params=p.Results;
    if ~isscalar(string(params.min_value_field))
        error('parameter "min_value_field" must be a character vector or string scalar.');
    end
    if ~isscalar(string(params.max_value_field))
        error('parameter "max_value_field" must be a character vector or string scalar.');
    end    
    nparam = height(T); 
    min_value = T.(params.min_value_field);    
    max_value = T.(params.max_value_field);
    if any(min_value>max_value)
        error('Values of "min_value" must be less than or equal to corresponding elements of "max_value"');
    end
    if ~istable(S)
        error('Input must be a cells table.');
    end
    ncells = height(S);
    %% apply inclusion criteria, looping over parameters
    include = true(ncells,1);
    table_vars=S.Properties.VariableNames;
    for i=1:nparam
        if ~ismember(T.parameter(i),table_vars)
            warning('%s is a unit inclusion parameter but is not a column in the cell table.',T.parameter(i));
            continue
        end
        out_of_range = S.(T.parameter(i)) < min_value(i) | S.(T.parameter(i)) > max_value(i);
        if any(out_of_range)
            include ( out_of_range ) = false;
            fprintf('%d units out of necessary range [%g,%g] for parameter %s.\n',sum(out_of_range),min_value(i),max_value(i),T.parameter(i));
        end
    end
    fprintf('%d of %d units (%.1f%%) DO NOT satisfy inclusion criteria.\n',sum(~include),numel(include),100*sum(~include)./numel(include));
end