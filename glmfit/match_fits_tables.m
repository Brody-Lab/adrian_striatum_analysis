function [table1,table2] = match_fits_tables(table1,table2)
    
    names{1} =  inputname(1);
    names{2} =  inputname(2);

    keys = {'cellno','recording_name'};
    if any(~ismember(keys,table1.Properties.VariableNames)) || any(~ismember(keys,table2.Properties.VariableNames))
        error('keys are not variable names in both tables.');
    end

    keyIdx1 = ismember(table1.Properties.VariableNames,keys);
    keyIdx2 = ismember(table2.Properties.VariableNames,keys);   
    
    l1 = ismember(table1(:,keyIdx1),table2(:,keyIdx2));
    l2 = ismember(table2(:,keyIdx1),table1(:,keyIdx2));
    
    if any(~l1)
        fprintf('Removed %d fits from table 1: %s.\n',sum(~l1),names{1});
    end
    
    if any(~l2)
        fprintf('Removed %d fits from table 2: %s.\n',sum(~l2),names{2});
    end     
    
    table1=table1(l1,:);
    table2=table2(l2,:);



end