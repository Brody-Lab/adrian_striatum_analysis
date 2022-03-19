%% load cells file for recording name  'A297_2021_07_18_02'
% Cells = load_Cells_file('A297_2021_07_18_02');
% 
% PB_set_constants;
% %% add spike_time_s field for clicks
% states2add  = {'left_clicks','right_clicks'};
% for c=1:numel(Cells.raw_spike_time_s)
%     for s=1:numel(states2add)
%         Cells.spike_time_s.(states2add{s}){c,1} = ...
%             group_spike_times(Cells.raw_spike_time_s{c,1}, Cells.Trials.stateTimes.(states2add{s}), kSpikeWindowS.(states2add{s}));
%     end
% end
% 
% %% exclude bad trials
% good_trials=find(~Cells.Trials.violated & ~isnan(Cells.Trials.stim_dur_s_actual));
% good_left_clicks = ismember(Cells.Trials.stateTimes.left_click_trial,good_trials);
% good_right_clicks = ismember(Cells.Trials.stateTimes.right_click_trial,good_trials);
% 
% %% make psths
% kPETH = get_PETH_params('std_s',0.01,'std_s_clicks',0.05,'type','LGAUSS');
% % these lines take about a minute each
% tic;psth_left=get_psth(Cells,1:numel(Cells.raw_spike_time_s),'states',{'left_clicks'},'kPETH',kPETH,'trial_idx',good_left_clicks,'nresamples',50);toc
% tic;psth_right=get_psth(Cells,1:numel(Cells.raw_spike_time_s),'states',{'right_clicks'},'kPETH',kPETH,'trial_idx',good_right_clicks,'nresamples',50);toc
% psth_right = psth_right(Cells.is_in_dorsal_striatum);
% psth_left = psth_left(Cells.is_in_dorsal_striatum);
%     

clear h;
%% plot example cells
id= kPETH.timeS.left_clicks<0.2 & kPETH.timeS.left_clicks>-0.01;
cellno = randsample([31  45 30 18 23  33 27  38],8);
%cellno = randsample(1:numel(psth_left),8);
figure('color',[1 1 1],'units','normalized','position',[0.4 0.3 0.26 0.65]);
for i=1:numel(cellno)
    b=subplot(4,2,i);
    set(b,P.axes_properties{:});
    set(gca,'xlim',[-0.03 0.35],'xgrid','on');     
    idx = kPETH.timeS.left_clicks> -0.01 & kPETH.timeS.left_clicks <= 0;
    mn = mean(psth_left(cellno(i)).left_clicks(:,idx));
    mn=mean(mn);
    line([-0.05 0.2],[1 1],'color',[1 1 1]/2,'LineWidth',1.5);
        line([ 0 0], [0.75 1.25],'color',[1 1 1]/2,'LineWidth',1.5);            
    h(1)=shadedErrorBar(kPETH.timeS.left_clicks(id),psth_left(cellno(i)).left_clicks(:,id)./mn,{@mean,@std},{'color',P.gamma_color_groups(end,:),'LineWidth',1.5});hold on;
    try
    h(1).patch.FaceAlpha=0.8;    
    end
    mn = mean(psth_right(cellno(i)).right_clicks(:,idx));
    mn=mean(mn);    
    h(2)=shadedErrorBar(kPETH.timeS.left_clicks(id),psth_right(cellno(i)).right_clicks(:,id)./mn,{@mean,@std},{'color',P.gamma_color_groups(2,:),'LineWidth',1.5});
    try
    h(2).patch.FaceAlpha=0.8;
    end
    if i==1
        legend([h.mainLine],{'left clicks','right clicks'});
    end
    %title((['A297\_2021\_07\_18\_02\_',num2str(cellno(i))]),'FontSize',8);
    if i==7
        xlabel('Time after click (s)');
        ylabel('Normalized Firing Rate');
    end
    %line(get(gca,'xlim'),[1 1],'color','k');
    %line([0 0],get(gca,'ylim'),'color','k','LineStyle','--');  
    set(gca,'xtick',0,'ytick',1);
    grid on;
    set(gca,'ylim',[0.5 1.5]);
    pos=get(gca,'position');
    if i/2 == round(i/2)
    pos(1) = pos(1)-0.12;
    else
        pos(1) = pos(1)+0.05;
    end
    if i>2
        pos(2) = pos(2)+0.076*ceil((i-2)/2);
    end
    set(gca,'position',pos,'FontSize',16,'LineWidth',1.3);
    if i==numel(cellno)
        line([0.27 0.27],[0.7 1.2],'color','k','LineWidth',1.5);
        line([0.17 0.27],[0.7 0.7],'color','k','LineWidth',1.5);
        text(0.2,0.6,'0.1 s','FontSize',13);
        text(0.34,0.95,{'50%','change'},'FontSize',15,'HorizontalAlignment','center');
        
    end
    %title(sprintf('%g',cellno(i)));
    axis off;    
    if i==numel(cellno)-1 
        text(-0.01,0.68,'0','FontSize',15);
        text(-0.05,1,'1','FontSize',15);
        text(-0.02,0.53,'Time after click','FontSize',15);
        text(-0.235,1,{'Normalized','Firing Rate'},'FontSize',15);
    else
        
    end

end