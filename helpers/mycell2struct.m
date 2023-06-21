function S = mycell2struct(C,mode)
    % C is a cell array where each element is a scalar struct
    
    % S is a struct array.
    
    % assuming fields of C are all identical, S(i) = C{i};
    % But in this case you could have just done S = [C{:}];
    
    % The advantage here is fields don't have to match.
    
    % Second argument tells what to do when they don't match.
    % 'union' means have the fields of S be the union, and have empties in
    % some places.
    % 'intersect' means throw out non-matching fields.
    
    fields={};
    if ~iscell(C)
        error('Input must be a cell array.');
    end
    if ndims(C)>2
        error('Input cell array cannot have more than two dimensions currently.');
    end
    validatestring(mode,{'union','intersect'});
    % pass once through C to validate elements and figure out fields of
    % output struct    
    for i=1:numel(C)
        if isscalar(C{i}) && isstruct(C{i})
            switch mode
                case 'union'
                    fields = union(fields,fieldnames(C{i}));
                case 'intersect'
                    fields = intersect(fields,fieldnames(C{i}));            
            end
        else
            error('All elements of input must be scalar structs.');
        end
    end
    
    % second pass through C to build up S
    for i=1:numel(C)
        for f=1:numel(fields)
            if isfield(C{i},fields{f})
                S(i,1).(fields{f}) = C{i}.(fields{f});
            end
        end
    end
    
    S = reshape(S,size(C));
   


end