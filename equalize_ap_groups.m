function ap_group = equalize_ap_groups(ap_group)
    for i=1:4
        idx{i} = find(ap_group==i);
        n(i) = numel(idx{i});
    end
    min_n = min(n);
    ap_group=zeros(size(ap_group));
    for i=1:4
        idx{i} = randsample(idx{i},min_n);
        ap_group(idx{i}) = i;
    end 
end