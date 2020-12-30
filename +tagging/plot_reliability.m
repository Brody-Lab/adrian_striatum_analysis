function plot_reliability(params)
    %% parse and validate inputs
    arguments
        params.min_dv (1,1) {mustBeNumeric,mustBePositive} = 2
        params.untagged_color (1,3) {mustBeNumeric,mustBeLessThanOrEqual(params.untagged_color,1),...
            mustBeGreaterThanOrEqual(params.untagged_color,0)} = [1 1 1]/2
        params.tagged_color (1,3) {mustBeNumeric,mustBeLessThanOrEqual(params.tagged_color,1),...
            mustBeGreaterThanOrEqual(params.tagged_color,0)} = spectrumRGB(473)
        params.nbins (1,1) {mustBeNumeric,mustBePositive} = 10      
        params.log_scale (1,1) logical = false
        params.max_days_implanted (1,1) = 300
    end
    P = get_parameters;
    %% load cells table
    T = readtable(P.cells_table_path);
    
    %% find indices of tagged and untagged cells
    tagged = T.reliability>0.9 & T.DV>params.min_dv & T.days_implanted<params.max_days_implanted & T.sess_date<datetime('2020-06-01'); % use reliability>0.9 as the metric for now
    untagged = T.reliability<0.9 & ~isnan(T.reliability) & T.DV>params.min_dv & T.days_implanted<params.max_days_implanted & T.sess_date<datetime('2020-06-01');     
    
    %% plot
    h(1)=histogram(T.reliability(untagged),'NumBins',params.nbins,'BinLimits',[0 1],'FaceColor',params.untagged_color);hold on;
    h(2)=histogram(T.reliability(tagged),'NumBins',params.nbins,'BinLimits',[0 1],'FaceColor',spectrumRGB(473));
    set(gca,P.axes_properties{:});
    if params.log_scale
        set(gca,'yscale','log','ytick',[1 10 100 1000],'yticklabel',[1 10 100 1000],'ylim',[8e-1 3e3]);      
    end    
    xlabel('Fraction of Trials with Firing Rate Increase');
    ylabel('Number of Cells');
    legend(h,{'Not Tagged','Tagged (Putative D2+ SPN)'},'Location','northwest');  
end