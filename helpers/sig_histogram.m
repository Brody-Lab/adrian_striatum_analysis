function [meancp,pval] = sig_histogram(cps,varargin)
    p=inputParser;
    p.KeepUnmatched=true;
    dividingLine = 0;
    p.addRequired('cps',@(x)validateattributes(x,{'numeric'},{}));
    p.addRequired('pvals',@(x)validateattributes(x,{'numeric'},{'nonnegative'}));
    p.addParameter('nbins',15,@(x)validateattributes(x,{'numeric'},{'positive','integer','scalar'}));
    p.parse(cps,varargin{:});
    [cps,pvals]=struct2var(p.Results,{'cps','pvals'});
    hold on
    if ~ismember('pvals',p.UsingDefaults)
        if any(size(pvals)~=size(cps)) 
            error('cps and pvals must be the same size.');
        end
        sig=histogram(cps(pvals<0.05));
        maxOffset=max(abs(minmax(cps)-dividingLine));
        sig.BinLimits=dividingLine + [-maxOffset maxOffset];
        sig.NumBins=nearestEven(p.Results.nbins);
        nonsig=histogram(cps(pvals>0.05));
        nonsig.BinLimits=dividingLine + [-maxOffset maxOffset];
        nonsig.NumBins=nearestEven(p.Results.nbins);
        xval=sig.BinEdges([1:end-1])+sig.BinWidth/2;
        sigval=sig.Values;
        nonsigval=nonsig.Values;
        xl=sig.BinLimits;
        cla;
        bar_handle=bar(xval',[sigval;nonsigval]','stacked');
        bar_handle(2).FaceColor=[1 1 1];
        bar_handle(1).FaceColor=[0 0 0];
    else
        h = histogram(cps);
        h.NumBins=p.Results.nbins;
        h.BinLimits=minmax(cps);
    end
    meancp=nanmean(cps);
    bootmean=bootstrp(1000,@nanmean,cps);
    pval=min(empirical_p(dividingLine,bootmean,'low'),empirical_p(dividingLine,bootmean,'high'));
    if pval<0.05
        sigstring='*';
        if pval<0.01
            sigstring='**';
        end
    else
        sigstring='n.s.';
    end
    line([dividingLine dividingLine],get(gca,'ylim'),'LineStyle',':','color',[1 1 1]/2);
    yl=get(gca,'ylim');
    set(gca,'xlim',xl);
    patch([-abs(diff(xl))/20 0 abs(diff(xl))/20]+meancp,[0 -abs(diff(yl))/20 0]+yl(2)-abs(diff(yl))/20,[1 1 1]/2);
    text(meancp,yl(2)-abs(diff(yl))/20,sigstring,'HorizontalAlignment','center','FontSize',15,'VerticalAlignment','baseline');
    title(['mean CP = ',num2str(round(meancp,3))]);
    box off;
    set(gca,'xgrid','on');
    xlabel('Choice Probability');
    ylabel('Number of Cells');    
end