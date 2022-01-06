%% scratch for plotting 1-on-1 meeting figs for carlos 1/6/22

% TS state space ('A297 - 2021-07-19 - 19076606841 - 5 PCs') by evidence
% strength, aligned to clicks on

psth_strong_left=get_pc_psth(stats,'trial_idx',stats.Trials.gamma<-2 & abs(stats.Trials.gamma)<50 ,'states',{'clicks_on'});
psth_weak_left=get_pc_psth(stats,'trial_idx',stats.Trials.gamma>-2 & stats.Trials.gamma<0  & abs(stats.Trials.gamma)<50  ,'states',{'clicks_on'});
psth_weak_right=get_pc_psth(stats,'trial_idx',stats.Trials.gamma<2 & stats.Trials.gamma>0 & abs(stats.Trials.gamma)<50  ,'states',{'clicks_on'});
psth_strong_right=get_pc_psth(stats,'trial_idx',stats.Trials.gamma>2 & abs(stats.Trials.gamma)<50 ,'states',{'clicks_on'});

baseline_psth=get_pc_psth(stats);


figure;
set(gca,'ColorOrder',copper(4));
data = cat(1,psth_strong_left,psth_weak_left,psth_weak_right,psth_strong_right);
for i=1:4
    h(i) = plot_pc_psth(data(i,:),'clicks_on','plotdim',3,'pcs',1:5,'time_lim_s',[-0.2 0.7]);hold on;
end
legend(h,{'strong left','weak left','weak right','strong right'})
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])


% left versus right clicks

%to equalize trials
% nleft = cellfun(@numel,stats.Trials.leftBups);
% nright = cellfun(@numel,stats.Trials.rightBups);
% nleft = nleft(~stats.Trials.violated);
% nright = nright(~stats.Trials.violated);
% wrong_trials = ~stats.Trials.is_hit;
% wrong_trials = wrong_trials(find(~stats.Trials.violated));
% 
% for i=1:20
%     equalized_hits = equalize_trials([nleft nright],wrong_trials,[5 5]);
%     correct_psth=get_pc_psth(stats,'trial_idx',stats.trial_idx(equalized_hits),'states',{'left_clicks','right_clicks'});
%     for k=1:numel(correct_psth)
%         left_clicks(i,k,:) = correct_psth(k).left_clicks;
%         right_clicks(i,k,:) = correct_psth(k).right_clicks;        
%     end
% end
% for k=1:numel(correct_psth)
%     correct_psth(k).left_clicks = squeeze(mean(left_clicks,1));
%     correct_psth(k).right_clicks = squeeze(mean(right_clicks,1));
% end




correct_psth=get_pc_psth(stats,'trial_idx',stats.Trials.is_hit,'states',{'left_clicks','right_clicks'});
wrong_psth=get_pc_psth(stats,'trial_idx',~stats.Trials.is_hit,'states',{'left_clicks','right_clicks'});



pcs=1:3;
d=3;
figure;
subplot(1,2,1);
h(1)=plot_pc_psth(correct_psth,'left_clicks','time_lim_s',[-0.025 0.2],'plotdim',d,'pcs',pcs);hold on
h(2)=plot_pc_psth(correct_psth,'right_clicks','time_lim_s',[-0.025 0.2],'plotdim',d,'pcs',pcs);
ax(1)=gca;
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
title('Correct Trials');
legend(h,{'contra clicks','ipsi clicks'});
subplot(1,2,2);
h(1)=plot_pc_psth(wrong_psth,'left_clicks','time_lim_s',[-0.025 0.2],'plotdim',d,'pcs',pcs);hold on
h(2)=plot_pc_psth(wrong_psth,'right_clicks','time_lim_s',[-0.025 0.2],'plotdim',d,'pcs',pcs);
ax(2)=gca;
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
title('Error Trials');
legend(h,{'contra clicks','ipsi clicks'});
linkaxes(ax);