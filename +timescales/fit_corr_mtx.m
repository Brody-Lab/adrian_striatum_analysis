function fit=fit_corr_mtx(corr_mtx,bin_size_ms,params)
    arguments
        corr_mtx 
        bin_size_ms (1,1) {isnumeric,mustBePositive,mustBeFinite}
    end
    if iscell(corr_mtx)
        corr_mtx = cat(3,corr_mtx{:});
        corr_mtx = nanmean(corr_mtx,3);
    end
    validateattributes(corr_mtx,{'numeric'},{'matrix','square'},'plot_corr_mtx','corr_mtx',1);
    autocorr = timescales.get_autocorr_from_corr_mtx(corr_mtx);
    xs = (1:numel(autocorr)) * bin_size_ms;
    fit = timescales.do_fit(xs,autocorr);    
end