function [h,fit]=plot_autocorr(autocorr,bin_size_ms,params)
    arguments
        autocorr
        bin_size_ms (1,1) {isnumeric,mustBePositive,mustBeFinite}
        params.fit (1,1) logical = true;
        params.color (1,3) double = [0 0 0];
    end        
    P = get_parameters;
    if isvector(autocorr)
        n_bins=numel(autocorr);
    else
        n_bins=size(autocorr,2); % X x N where X is number of cells, N is number of bins
    end
    xs = (1:n_bins) * bin_size_ms;
    if isvector(autocorr)
        h=plot(xs,diags,'o:','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','none');
    else
        boots = bootstrp(1000,@nanmean,autocorr);
        autocorr=nanmean(autocorr);
        h=shadedErrorBar(xs,boots,{@mean,@std});        
        h.mainLine.MarkerFaceColor = params.color;
        h.mainLine.MarkerEdgeColor='none';
        h.mainLine.Marker = 'o';
        h.mainLine.LineStyle = ':';
        h.mainLine.MarkerSize=10;
        h.patch.FaceColor = params.color;
        h.patch.FaceAlpha=0.2;
    end
    set(gca,P.axes_properties{:});
    xlabel('lag (ms)');
    ylabel('Autocorrelation');
    if params.fit
        fit = timescales.do_fit(xs,autocorr);
        hold on;
        plot(xs,fit.f(fit.b,xs),'color',params.color);
    end
end