function is_in_dorsal_striatum = get_dorsal_striatal_cells(T)
    striatum_names={'TS','DMS','DLS','ADS','Str2','Str3'};
    if istable(T) % cells table
        is_in_dorsal_striatum = false(height(T),1);
        for i=1:height(T)
            if ~isempty(T.regions{i}) && contains(T.region{i},striatum_names)
                is_in_dorsal_striatum(i)=true;
            end
        end
    elseif isstruct(T) % cells file (already processed one in db that has up to date penetrations fields)
        ncells = numel(T.region);
        is_in_dorsal_striatum = is_in_region(T,1:ncells,striatum_names);% calls new function "is_in_region" which encapsulates useful code and works automatically on joined files - AGB, 8/9/2022
    end
end