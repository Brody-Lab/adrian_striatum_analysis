function plot_latency(params)
    %% parse and validate inputs
    arguments
        params.min_dv (1,1) {mustBeNumeric,mustBePositive} = 2
        params.nan_offset_s (1,1) {mustBeNumeric,mustBePositive} = 0.015
        params.log_scale (1,1) logical = true
        params.untagged_color (1,3) {mustBeNumeric,mustBeLessThanOrEqual(params.untagged_color,1),...
            mustBeGreaterThanOrEqual(params.untagged_color,0)} = [1 1 1]/2
        params.tagged_color (1,3) {mustBeNumeric,mustBeLessThanOrEqual(params.tagged_color,1),...
            mustBeGreaterThanOrEqual(params.tagged_color,0)} = spectrumRGB(473)
        params.nbins (1,1) {mustBeNumeric,mustBePositive} = 35        
    end
    P = get_parameters;
    %% load cells table
    T = readtable(P.cells_table_path);
    
    %% find indices of tagged and untagged cells
    tagged = T.reliability>0.9; % use reliability>0.9 as the metric for now
    untagged = T.reliability<0.9 & ~isnan(T.reliability);    
    tagged_times = T.first_sig_time_s(tagged & T.DV>params.min_dv);    
    untagged_times = T.first_sig_time_s(untagged & T.DV>params.min_dv);  
    
    %% convert NaNs to special value for plotting
    max_time = max(max(tagged_times),max(untagged_times));
    tagged_times(isnan(tagged_times))= max_time+params.nan_offset_s;
    untagged_times(isnan(untagged_times))= max_time+params.nan_offset_s;    
    
    %% plot
    h(1)=histogram(untagged_times,'NumBins',params.nbins,'BinLimits',[0 max_time+params.nan_offset_s],'FaceColor',params.untagged_color);hold on;
    h(2)=histogram(tagged_times,'NumBins',params.nbins,'BinLimits',[0 max_time+params.nan_offset_s],'FaceColor',spectrumRGB(473));
    set(gca,P.axes_properties{:},'xtick',[0:0.01:max(max_time,0.1) max_time+params.nan_offset_s],...
        'xticklabel',[arrayfun(@num2str,1000*(0:0.01:max(max_time,0.1)),'uni',0) 'N/A']);
    if params.log_scale
        set(gca,'yscale','log','ytick',[1 10 100 1000],'yticklabel',[1 10 100 1000],'ylim',[8e-1 3e3]);      
    end
    xlabel('Time (ms) to Significant Firing Rate Increase');
    ylabel('Number of Cells');
    legend(h,{'Not Tagged','Tagged (Putative D2+ SPN)'},'Location','northwest');  
end