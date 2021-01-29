function fit=fit_corr_mtx(corr_mtx,bin_size_ms)
    arguments
        corr_mtx 
        bin_size_ms (1,1) {isnumeric,mustBePositive,mustBeFinite}
    end
    if iscell(corr_mtx)
        corr_mtx = cat(3,corr_mtx{:});
        corr_mtx = nanmean(corr_mtx,3);
    end
    validateattributes(corr_mtx,{'numeric'},{'matrix','square'},'plot_corr_mtx','corr_mtx',1);
    n=size(corr_mtx,1);  
    autocorr = get_autocorr_from_corr_mtx(corr_mtx);
    diags= spdiags(corr_mtx);
    nans = triu(spdiags(NaN(n)));
    diags(~isnan(nans))=NaN;
    diags=nanmean(diags(:,n+1:end));
    xs = (1:numel(diags)) * bin_size_ms;
    fit = timescales.do_fit(xs,diags);    
end