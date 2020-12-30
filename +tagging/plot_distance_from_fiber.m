function plot_distance_from_fiber(params)
    %% parse and validate inputs
    arguments
        params.min_dv (1,1) {mustBeNumeric,mustBePositive} = 1.5
        params.nbins (1,1) {mustBeNumeric,mustBePositive} = 4
    end
    P = get_parameters;
    %% load cells table
    T = readtable(P.cells_table_path);
    
    %% find indices of tagged and untagged cells
    tagged = T.reliability>0.9 & T.DV>params.min_dv; % use reliability>0.9 as the metric for now
    untagged = T.reliability<0.9 & ~isnan(T.reliability) & T.DV>params.min_dv;    
    distance_from_fiber = T.distance_from_fiber_tip;
    
    %% count bins
    bin_limits = [min(distance_from_fiber(tagged|untagged)) max(distance_from_fiber(tagged|untagged))];
    bin_edges = linspace(bin_limits(1),bin_limits(2),params.nbins+1);
    for i=params.nbins:-1:1
        in_range = distance_from_fiber>bin_edges(i) & distance_from_fiber<bin_edges(i+1);
        n_tagged(i) = sum(tagged & in_range) ;
        total_recorded(i) = sum((tagged|untagged) & in_range);
        bin_mids(i) = mean(bin_edges([i i+1]));
    end
    
    %% plot
    [~,pci] = binofit(n_tagged,total_recorded,0.05);
    frac_tagged = n_tagged./total_recorded;
    bar(bin_mids,100*frac_tagged);hold on;
    h=errorbar(bin_mids,100*frac_tagged,100*abs(frac_tagged-pci(:,1)'),100*abs(frac_tagged-pci(:,2)'));
    h.Color=[0 0 0];h.LineStyle='none';
    xl=xlim;
    set(gca,P.axes_properties{:},'xlim',[0 xl(2)]);
    xlabel('Distance from Fiber Tip (mm)');
    ylabel('Percent Tagged');
end