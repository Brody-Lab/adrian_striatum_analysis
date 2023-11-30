function plot_pref_click_amplitude(stats,pref_mode,varargin)
    % uses max abs deviation to determine click preference. and the sign of
    % the prefered click's max abs deviation (i.e. < or >1) determines
    % whether it's treated as a negative click of not).
    
    % fill group_var in with NaNs if you don't want to include some
    % entries. i.e. some entries are not in any group.
    
    %% parse and validate inputs
    p=inputParser;
    P=get_parameters();
    validateattributes(stats,{'struct'},{'nonempty','vector'},'plot_average_click_by_preference','stats',1);    
    p.addRequired('pref_mode',@(x)validateattributes(x,{'string'},{'nonempty','scalar'}));
    p.addParameter('group_by',[],@(x)validateattributes(x,{'numeric','logical'},{'vector'}));
    p.addParameter('separate_negative',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('flip_negative',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('equalize_groups',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('plot_average',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.parse(pref_mode,varargin{:});
    params=p.Results;
    group_var = params.group_by;
    stats=stats(:);
    N= numel(stats);
    if ismember('group_by',p.UsingDefaults)
        group_var = ones(N,1);
    elseif numel(group_var)~=N
        error('length of group_var must match number of elements of stats.');
    else
        group_var=group_var(:);
    end
    [group_var,sort_idx] = sort(group_var);
    stats=stats(sort_idx);
    validatestring(pref_mode,["deviation","average","both"],'plot_average_click_by_preference','pref_mode',2);
    [g_idx,id] = findgroups(group_var);
    if params.separate_negative && params.flip_negative
        warning('flip negative has no effect when separate negative is true.');
    end
        
    %% extract relevant data from stats
    left_clicks = get_combined_weights(stats,'left_clicks');
    right_clicks = get_combined_weights(stats,'right_clicks');
    covariate_stats = {stats.covariate_stats};
    
    avg_clicks = right_clicks;
    avg_clicks.data = (right_clicks.data + left_clicks.data) /2;    
    switch pref_mode 
        case "deviation"
            is_negative = cellfun(@(x)x.pref_click_max_deviation,covariate_stats)<1;            
        case "average" 
            is_negative = cellfun(@(x)x.pref_click_average,covariate_stats)<1;     
        case "both"
            is_negative = cellfun(@(x)x.pref_click_average,covariate_stats)<1;  % arbitrary choice                  
    end  
    tr=buildGLM.get_tr(avg_clicks.tr);
    colors=P.ap_group_colors; 
    for i=1:numel(id)
        n(i)=sum(g_idx==i);
    end
    for i=1:numel(id)
        if params.equalize_groups && numel(id)>1 && n(i)>min(n)
            idx{i} = randsample(find(g_idx==i),min(n));
        else
            idx{i} = find(g_idx==i);
        end
    end 
    set(gca,'xlim',[0 1],'box','off');
    if params.separate_negative
        subplot(2,1,1);
        for i=1:numel(id)
            h(i) = shadedErrorBar(tr, (exp((avg_clicks.data(intersect(idx{i},find(is_negative)),:)))),{@mean,@(x)std(x)./sqrt(size(x,1))},{'color',colors(i,:),'LineWidth',1.5});hold on;
            h(i).patch.FaceAlpha=0.5;            
        end
        subplot(2,1,2);        
        for i=1:numel(id)
            h(i)=shadedErrorBar(tr, (exp((avg_clicks.data(intersect(idx{i},find(~is_negative)),:)))),{@mean,@(x)std(x)./sqrt(size(x,1))},{'color',colors(i,:),'LineWidth',1.5});hold on;
            h(i).patch.FaceAlpha=0.5;
        end    
        matchylim(gcf);    
    else
       avg_clicks.data(is_negative==1,:) =   -avg_clicks.data(is_negative==1,:);
       if params.plot_average
            for i=1:numel(id)
                h(i) = shadedErrorBar(tr, ((avg_clicks.data( idx{i},:))),{@(x)exp(mean(x)),@(x)std(x)./sqrt(size(x,1))},{'color',colors(i,:),'LineWidth',1.5});hold on;  
                h(i).patch.FaceAlpha=0.5;            
            end    
       else
           if params.equalize_groups
               N = min(n)*numel(id);
           end
          figure;imagesc(1:numel(tr),1:N,exp(avg_clicks.data(unique(cat(1,idx{:})),:)));set(gca,'clim',[0.7 1.4]);set(gca,'xscale','linear');
       end
        
        %matchylim(gcf);        
    end
    l=legend([h.mainLine],P.ap_group_labels);
    l.Box='off';
    ylabel('Gain');
    xlabel('Time after click (s)');
    set(gca,'ytick',[1 1.05 1.1],P.axes_properties{:});
end




