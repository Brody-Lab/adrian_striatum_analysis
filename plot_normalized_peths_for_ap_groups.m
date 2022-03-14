clear psths_norm
for i=1:numel(P.ap_groups)
    is_in_group = cells_table.AP>=P.ap_groups{i}(1) & cells_table.AP<P.ap_groups{i}(2);
    cells_table.ap_group(is_in_group) = i;
    
    psths{i} = cat(3,cells_table.psth_clicks_on{cells_table.ap_group==i});
    pref = cells_table.pref_left_MI(cells_table.ap_group==i);
    pref_pval = cells_table.pref_left_MI_pval(cells_table.ap_group==i);    
    
    no_pref = isnan(pref) | isnan(pref_pval) | pref_pval>0.05;
    
    clear bad rate flips
    for k=1:size(psths{i},3)
        this_one = psths{i}(:,1:100,k);
        rate(k) = nanmean(this_one(:));
        if rate(k)<0.05
            bad(k)=true;
        else
            bad(k)=false;
        end
        psths_norm{i}(:,:,k) = psths{i}(:,:,k) ./ rate(k);

    end
    
    psths_norm{i} = psths_norm{i}(:,:,~bad' & ~no_pref );
    pref= pref(~bad' & ~no_pref);
    psths_norm{i}(:,:,pref>0) = flip(psths_norm{i}(:,:,pref>0));
    
end



figure('units','normalized','position',[0.07 0 0.19 1.1],'color',[1 1 1]);
for i=1:4
    h(i) = subplot(4,1,i);set(h(i),'units','normalized');
    idx = kPETH.timeS.clicks_on==0;
    this_psth = bsxfun(@minus,psths_norm{i},psths_norm{i}(:,idx,:));
    %this_psth = psths_norm{i}-1;
    this_psth = this_psth(:,20:10:end,:);
    err = nanstd(this_psth,[],3) ./ sqrt(size(this_psth,3));
    %line([-1 1],[0 0],'color',[0 0 0 ]+0.5,'linewidth',2);hold on;
    for k=1:4
        g(k)=errorbar(kPETH.timeS.clicks_on(20:10:end),nanmean(this_psth(k,:,:),3),err(k,:),'LineWidth',1.5);colororder(P.gamma_color_groups);hold on;
        g(k).CapSize=0;
    end
    gammas={'strong preferred','weak preferred','weak null','strong null'};
    if i==1
        for k=1:4
            text(0.02,2.6-0.3*k,gammas{k},'color',P.gamma_color_groups(k,:),'FontSize',14);
        end
    end
    set(gca,'xlim',[0 0.8],P.axes_properties{:});
    yl{i}=get(gca,'ylim');
    set(gca,'ylim',[-0.5 min(max(yl{i}(2),0.5),2.5)],...
        'ytick',[0.5 1 1.5 2 2.5 3 3.5]-1,'xtick',[0:0.2:1],'XGrid','on','YGrid','on','yticklabel',[0.5 1 1.5 2 2.5 3 3.5],'xticklabel',[0 0.2 0.4 0.6 0.8]);
    if i==4
        l=xlabel('Time after first click (s)');
        l.Position = [0.65 -0.9 -1];
    else
        set(gca,'xticklabel',[]);
    end
    if i==4
        ylabel({'Normalized','Firing Rate'});
    end
end

for i=1:4
    subplot(4,1,i)    ;drawnow;
    if i==1
        %text(0.02,1.5,{'Normalized','Firing','Rate'},'HorizontalAlignment','left','FontSize',14);
    end
    title(P.ap_group_labels{i},'color',P.ap_group_colors(i,:));
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
        pos(2) = pos(2)+0.065;
    end
    set(h(i),'units','normalized',P.axes_properties{:},'position',pos,'XAxisLocation','origin');
    
    
        
end
