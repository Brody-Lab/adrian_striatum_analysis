%% generated PETH plots for GRS talk
%% edited/improved 11/2022 before SfN

P=get_parameters();
ref_event = 'first_click';
mask_states = {'' 'cpoke_req_end'};
column_name = 'clicks_on_peth';
a=tic;
[cells_table,kPETH] = add_peth_to_cells_table(load_cells_table(),'ref_event',ref_event,...
    'separate_by','signal_strength','mask_states',mask_states,'column_name',column_name,'mask_state_offset',[0 0]);
fprintf('total loop took %s',timestr(toc(a)));
baseline_idx = kPETH.timeS.(ref_event)<0;
zrpt_idx = kPETH.timeS.(ref_event)==0;
accum_idx = kPETH.timeS.(ref_event)>0 & kPETH.timeS.(ref_event)<0.5;


% loop over cells to find baseline rate by which to normalize
ncells = height(cells_table);
baseline = zeros(ncells,1);
accum_rate = zeros(ncells,1);
psths = cell(ncells,1);
for i=1:ncells
    tmp = cells_table.(column_name){i}(:,accum_idx);
    accum_rate(i) = nanmean(tmp(:));    
    tmp = cells_table.(column_name){i}(:,baseline_idx);
    baseline(i) = nanmean(tmp(:));    
    psths{i} = cells_table.(column_name){i} ./ baseline(i); % divide by baseline
    psths{i} = psths{i}(:,20:10:end) - psths{i}(:,zrpt_idx) +1 ; % normalize so that values are on average 1 at time 0
    if cells_table.choice_MI(i)>0
        psths{i} = flip(psths{i}); % flip left-preferring so that psths can be interpreted as pref/null
    end
end

figure('units','normalized','color',[1 1 1]);
for i=1:numel(P.ap_groups)
    is_in_group = cells_table.ap_group==i & cells_table.is_in_dorsal_striatum;
    sig_pref_idx = is_in_group & ~isnan(cells_table.choice_MI) & cells_table.choice_MI_pval<0.05 & baseline>0.05 ; % find indices for significantly choice-preferring neurons
    sum(sig_pref_idx)
    psth_data = cat(3,psths{sig_pref_idx});
    h(i) = subplot(1,4,i);set(h(i),'units','normalized');
    for k=1:(numel(P.gamma_ranges)-1)
        err = bootstrp(1000,@nanmean,squeeze(psth_data(k,:,:))');        
        g(k)=errorbar(kPETH.timeS.(ref_event)(20:10:end),nanmean(psth_data(k,:,:),3),std(err) ,'LineWidth',2);
        colororder(P.gamma_color_groups);hold on;
        g(k).CapSize=0;
    end    
    if i==1
        for k=1:4
            text(0.02,3.5-0.2*k,P.gamma_labels{k},'color',P.gamma_color_groups(k,:),'FontSize',14);
        end
    end
    set(gca,'xlim',[0 0.8],P.axes_properties{:});
    set(gca,'ylim',[0.8 3.5],...
        'ytick',0:5,'xtick',[0:0.2:1],'XGrid','off',...
        'YGrid','off','xticklabel',[0 0.2 0.4 0.6 0.8]);
        l=xlabel('Time after first click (s)');

    if i==1
        ylabel({'Normalized','Firing Rate'});
    end    
end




for i=1:4
    subplot(4,1,i)    ;drawnow;
    if i==1
        %text(0.02,1.5,{'Normalized','Firing','Rate'},'HorizontalAlignment','left','FontSize',14);
    end
    d{i} = daspect;
    pb{i} = pbaspect;    
    yl{i} = get(gca,'ylim');    
    yratio = range(yl{i})./range(yl{1});
    pb{i}(2) = pb{i}(2).*yratio;
    pbaspect(pb{i});
    %daspect(d{i});    
    h(i) = subplot(4,1,i);
    pos=h(i).Position;
    pos(2) = pos(2)+0.08*(i-1);
    if i==3
        pos(2)=pos(2)+0.04;
    elseif i==4
        pos(2) = pos(2)+0.055;
    end
    pos(1)=pos(1)+0.05;
    set(h(i),'units','normalized',P.axes_properties{:},'position',pos,'XAxisLocation','origin','FontSize',16,'LineWidth',1.4);
    title(P.ap_group_labels{i},'color',P.ap_group_colors(i,:),'FontSize',21);    
    grid off;
    
    
        
end
