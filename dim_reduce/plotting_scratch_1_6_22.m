%% scratch for plotting 1-on-1 meeting figs for carlos 1/6/22

% TS state space ('A297 - 2021-07-19 - 19076606841 - 5 PCs') by evidence
% strength, aligned to clicks on


psth_strong_left=get_pc_psth(stats,'trial_idx',stats.Trials.gamma<-2 & abs(stats.Trials.gamma)<50 & stats.Trials.is_hit ,'states',{'clicks_on','cpoke_out'});
psth_weak_left=get_pc_psth(stats,'trial_idx',stats.Trials.gamma>-2 & stats.Trials.gamma<0  & abs(stats.Trials.gamma)<50 & stats.Trials.is_hit ,'states',{'clicks_on','cpoke_out'});
psth_weak_right=get_pc_psth(stats,'trial_idx',stats.Trials.gamma<2 & stats.Trials.gamma>0 & abs(stats.Trials.gamma)<50 & stats.Trials.is_hit ,'states',{'clicks_on','cpoke_out'});
psth_strong_right=get_pc_psth(stats,'trial_idx',stats.Trials.gamma>2 & abs(stats.Trials.gamma)<50 & stats.Trials.is_hit ,'states',{'clicks_on','cpoke_out'});

baseline_psth=get_pc_psth(stats,'trial_idx',abs(stats.Trials.gamma)<50 & stats.Trials.is_hit & z,'states',{'clicks_on','cpoke_out'});



pcs = [2 3 5]; % A297 FOF or ADS
pcs = [1 2 3]; %A294 FOF
pcs = [1 2 4]; %A294 ADS
pcs = [ 1 3 4]; % A297 TS (all regions)

figure;
set(gca,'ColorOrder',copper(4));
data = cat(1,psth_strong_left,psth_weak_left,psth_weak_right,psth_strong_right);
for i=1:4
    h(i) = plot_pc_psth(data(i,:),'clicks_on','plotdim',3,'pcs',pcs,'time_lim_s',[0 0.75],'start_at_zero',true,'baseline',baseline_psth);hold on;
end
legend(h,{'strong left','weak left','weak right','strong right'})
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])

% cpoke out

figure;
set(gca,'ColorOrder',copper(4));
data = cat(1,psth_strong_left,psth_weak_left,psth_weak_right,psth_strong_right);
for i=1:4
    h(i) = plot_pc_psth(data(i,:),'cpoke_out','plotdim',3,'pcs',[1 2 4],'time_lim_s',[-1.5 1.5]);hold on;
end
legend(h,{'strong left','weak left','weak right','strong right'})
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% set(gca,'ztick',[])


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



pcs=[1 3 4]; % A297
pcs = [2 4 5]; %A242
pcs=[1 2 5];
d=3;
figure;
subplot(1,2,1);
h(1)=plot_pc_psth(correct_psth,'left_clicks','time_lim_s',[-0.025 0.2],'plotdim',d,'pcs',pcs,'start_at_zero',true);hold on
h(2)=plot_pc_psth(correct_psth,'right_clicks','time_lim_s',[-0.025 0.2],'plotdim',d,'pcs',pcs,'start_at_zero',true);
ax(1)=gca;
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
title('Correct Trials');
legend(h,{'contra clicks','ipsi clicks'});
subplot(1,2,2);
h(1)=plot_pc_psth(wrong_psth,'left_clicks','time_lim_s',[-0.025 0.2],'plotdim',d,'pcs',pcs,'start_at_zero',true);hold on
h(2)=plot_pc_psth(wrong_psth,'right_clicks','time_lim_s',[-0.025 0.2],'plotdim',d,'pcs',pcs,'start_at_zero',true);
ax(2)=gca;
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
title('Error Trials');
legend(h,{'contra clicks','ipsi clicks'});
linkprop(ax,{'xlim','ylim','zlim'});