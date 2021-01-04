function plot_metrics(params)
    %% parse and validate inputs
    arguments
        params.min_dv (1,1) {mustBeNumeric,mustBePositive} = 1.5
        params.untagged_color (1,3) {mustBeNumeric,mustBeLessThanOrEqual(params.untagged_color,1),...
            mustBeGreaterThanOrEqual(params.untagged_color,0)} = [1 1 1]/2
        params.tagged_color (1,3) {mustBeNumeric,mustBeLessThanOrEqual(params.tagged_color,1),...
            mustBeGreaterThanOrEqual(params.tagged_color,0)} = spectrumRGB(473)
        params.nbins (1,1) {mustBeNumeric,mustBePositive} = 30      
        params.log_scale (1,1) logical = false
    end
    P = get_parameters;
    %% load cells table
    T = readtable(P.cells_table_path);
    
    include = T.DV>params.min_dv;
    
    %% find indices of tagged and untagged cells
    tagged = T.reliability>0.9 & T.DV>params.min_dv; % use reliability>0.9 as the metric for now
    untagged = T.reliability<0.9 & ~isnan(T.reliability) & T.DV>params.min_dv ;     
    
    fixfun = @(x)x;
    
    %% plot
    count=0;
    k=0;
    l=0;
    n=length(P.tagging_metrics);
    for i=1:n
        for j=1:n
            count=count+1;
            if i==j
                subplot(n,n,count);          
                histogram(fixfun(T.(P.tagging_metrics{i})(include)),'NumBins',params.nbins);
            elseif i>j
                subplot(n,n,count);               
                scatter(fixfun(T.(P.tagging_metrics{j})(include)),fixfun(T.(P.tagging_metrics{i})(include)));                
            end
            if mod(count,n)==1
                k=k+1;
                h=ylabel(P.tagging_metric_names{k});h.Rotation=0;h.HorizontalAlignment='right'
            end
            if count>(n-1)*n
                l=l+1;
                xlabel(P.tagging_metric_names{l});
            end                
        end
    end
end