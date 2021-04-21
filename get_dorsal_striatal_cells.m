function is_in_dorsal_striatum = get_dorsal_striatal_cells(T)
    if istable(T) % cells table
        is_in_dorsal_striatum = false(height(T),1);
        for i=1:height(T)
            if ~isempty(T.regions{i}) && any(ismember({'TS','DMS','DLS'},T.regions{i}))
                is_in_dorsal_striatum(i)=true;
            end
        end
    elseif isstruct(T) % cells file
        T = import_penetration(T);
        ncells = numel(T.regions);
        is_in_dorsal_striatum = false(ncells,1);
        for i=1:ncells
            if T.regions(i)
                names = T.penetration.regions(T.regions(i)).name;
                if any(ismember({'TS','DMS','DLS'},names))
                    is_in_dorsal_striatum(i)=true;
                end
            end
        end
    end
end