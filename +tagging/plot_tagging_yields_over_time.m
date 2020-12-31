function plot_tagging_yields_over_time(params)
    %% parse and validate inputs
    arguments
        params.min_dv (1,1) {mustBeNumeric,mustBePositive} = 1.5
        params.max_fiber_distance (1,1) {mustBeNumeric,mustBePositive} = 5
    end
    P = get_parameters;
    %% load cells table
    T = readtable(P.cells_table_path);
    
    %% find indices of tagged and untagged cells
    tagged = T.reliability>0.9 & T.DV>params.min_dv & T.distance_from_fiber_tip<params.max_fiber_distance; % use reliability>0.9 as the metric for now
    untagged = T.reliability<0.9 & ~isnan(T.reliability) & T.DV>params.min_dv & T.distance_from_fiber_tip<params.max_fiber_distance;    
    
    [unique_sessids,first_idx,unique_idx] = unique(T.sessid);

    
    for i=1:length(unique_sessids)
        days_since_viral_injection(i,1) = T.days_since_viral_injection(first_idx(i));
        power(i,1) = T.laser_power_mW(first_idx(i));
        n_tagged(i,1) = sum(tagged(unique_idx==i)) ;
        n_recorded(i,1) = sum(tagged(unique_idx==i) | untagged(unique_idx==i));
        frac_tagged(i,1) = n_tagged(i)./n_recorded(i);
        rat(i,1) = T.rat(first_idx(i));
    end
    exclude_idx = isnan(power);
    
    frac_tagged(exclude_idx)=[];
    n_recorded(exclude_idx)=[];
    n_tagged(exclude_idx)=[];
    power(exclude_idx)=[];
    days_since_viral_injection(exclude_idx)=[];
    rat(exclude_idx)=[];
    
    group_vars = {'rat','power'};
    [u,uu,group_idx] = unique(table(rat,power));
    for i=1:height(u)
        these_idx = group_idx==i;
        if u.power(i)==1
            s=0.2;
        elseif u.power(i) == 5
            s=0.52;
        elseif u.power(i) == 10
            s=1;
        end
        if u.rat(i)=="A249"
            hue=0.55;
        else
            hue=0;
        end
        color=hsv2rgb([hue s 1]);
        h(i)=plot(days_since_viral_injection(these_idx),100*frac_tagged(these_idx),'o:','Color',color,'MarkerFaceColor',color);hold on
        labels{i} = sprintf('Rat %s, %g mW',rat{uu(i)},power(uu(i)));
    end
    xl=xlim;
    set(gca,P.axes_properties{:},'xlim',[0 xl(2)]);
    yl=ylim;    
    xlabel('Days Since Viral Injection');
    ylabel('Percent Tagged Cells');
    patch([65 245 245 65],yl([1 1 2 2]),[ 0 0 0],'FaceAlpha',0.1,'EdgeColor','none');
    legend(h,labels,'Location','southwest');    
    set(gca,'ylim',yl);
    text(110,25,{'COVID-19 Disruption'});
end