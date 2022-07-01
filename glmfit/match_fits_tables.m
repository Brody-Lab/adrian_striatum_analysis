function [l1,l2] = match_fits_tables(table1,table2)
    
    keys = {'cellno','recording_name'};
    if any(~ismember(keys,table1.Properties.VariableNames)) || any(~ismember(keys,table2.Properties.VariableNames))
        error('keys are not variable names in both tables.');
    end

    keyIdx1 = ismember(table1.Properties.VariableNames,keys);
    keyIdx2 = ismember(table2.Properties.VariableNames,keys);
    
    l1 = ismember(table1(:,keyIdx1),table2(:,keyIdx2));
    l2 = ismember(table2(:,keyIdx1),table1(:,keyIdx2));



end