function [fits_table,maxi] = get_best_fit(fits_tables)
    validateattributes(fits_tables,{'cell'},{'vector','nonempty'},mfilename,'fits_tables',1);
        
    for i=1:numel(fits_tables)
        if i<numel(fits_tables)
            [fits_tables{i},fits_tables{i+1}] = match_fits_tables(fits_tables{i},fits_tables{i+1});
        end
        if i>1
            good = good | [fits_tables{i}.stats.success];
        else
            good = [fits_tables{i}.stats.success];
        end
    end
    
    for i=numel(fits_tables):-1:1
        fits_tables{i} = fits_tables{i}(good,:);
        metrics{i} = [fits_tables{i}.stats.metrics];
        bps(:,i) = [metrics{i}.bits_per_spk];
    end
    
    [~,maxi] = max(bps,[],2);
    
    fits_table = fits_tables{1};
    for i=2:numel(fits_tables)
       idx = maxi==i;
       fits_table(idx,:) = fits_tables{i}(idx,:);
    end
end