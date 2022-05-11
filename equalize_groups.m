function ap_group = equalize_groups(ap_group,min_n)
    u=unique(ap_group);
    nu=numel(u);
    for i=1:nu
        idx{i} = find(ap_group==u(i));
        n(i) = numel(idx{i});
    end
    if nargin<2
        min_n = min(n);
    end 
    ap_group=NaN(size(ap_group));
    for i=1:nu
        idx{i} = randsample(idx{i},min_n);
        ap_group(idx{i}) = i;
    end 
end