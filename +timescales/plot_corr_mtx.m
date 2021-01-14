function [h,fit]=plot_corr_mtx(corr_mtx,bin_size_ms,params)
    arguments
        corr_mtx 
        bin_size_ms (1,1) {isnumeric,mustBePositive,mustBeFinite}
        params.mode string {ismember(params.mode,{'full','diag'})} = 'diag';
        params.fit (1,1) logical = true;
    end
    P = get_parameters;
    if iscell(corr_mtx)
        corr_mtx = cat(3,corr_mtx{:});
        corr_mtx = nanmean(corr_mtx,3);
    end
    validateattributes(corr_mtx,{'numeric'},{'matrix','square'},'plot_corr_mtx','corr_mtx',1);
    n=size(corr_mtx,1);  
    switch params.mode
        case 'diag'
            diags= spdiags(corr_mtx);
            nans = triu(spdiags(NaN(n)));
            diags(~isnan(nans))=NaN;
            diags=nanmean(diags(:,n+1:end));
            xs = (1:numel(diags)) * bin_size_ms;
            h=plot(xs,diags,'o:','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','none');
            set(gca,P.axes_properties{:});
            xlabel('lag (ms)');
            ylabel('Autocorrelation');
            if params.fit
                fit = timescales.do_fit(xs,diags);
                hold on;
                plot(xs,fit.f(fit.b,xs));
            end
        case 'full'
    end
end