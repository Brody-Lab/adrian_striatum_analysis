function autocorr = get_autocorr_from_corr_mtx(corr_mtx)
    n=size(corr_mtx,1);
    autocorr= spdiags(corr_mtx,1:(n-1)); 
    nans = spdiags(NaN(n),1:(n-1));
    autocorr(~isnan(nans))=NaN;
    autocorr=nanmean(autocorr);
end