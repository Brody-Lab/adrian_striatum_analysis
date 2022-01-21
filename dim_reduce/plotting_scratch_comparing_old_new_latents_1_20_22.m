file_old=load('C:\Users\abondy\Downloads\Bondy_Post_Striatum.mat'); % 0-3.25 slow latents
file=load('C:\Users\abondy\Downloads\Bondy_Post_Striatum_fast.mat'); %0-1.5 seconds fast latents

Cells_TS = load('X:\RATTER\PhysData\NP_sorted\Adrian\A297\A297_2021_07_19_02\A297_2021_07_19_02_g0\spikesort_2021_09_21_16_19_19_ks2jrc\A297_2021_07_19_02_g0_t0.imec0.ap_res.Cells.mat');

gpfa_old=format_latents(Cells_TS,file_old.latents,'cpoke_in',[0 3.25],'trial_idx',~Cells_TS.Trials.violated);
gpfa=format_latents(Cells_TS,file.latents,'cpoke_in',[0 1.5],'trial_idx',~Cells_TS.Trials.violated);

kPETH = get_PETH_params('type','LGAUSS','std_s_clicks',0.015,'std_s',0.015,'resolution_s',5e-3,'resolution_s_clicks',5e-3);


[gpfa_old.psth,gpfa_old.psth_std] = get_pc_psth(gpfa_old,'nresamples',25,'states',{'clicks_on'},'kPETH',kPETH);
[gpfa.psth,gpfa.psth_std] = get_pc_psth(gpfa,'nresamples',25,'states',{'clicks_on'},'kPETH',kPETH);


figure;
for i=1:5
    subplot(5,1,i);
    shadedErrorBar(kPETH.timeS.clicks_on,gpfa.psth(i).clicks_on,gpfa.psth_std(i).clicks_on );hold on;
    xlabel('Time (s) since stimulus start');
    ylabel(sprintf('Factor %g',i));
    yl=get(gca,'ylim');
    h=line([0 0],yl);
    h.LineStyle=':';
    h.LineWidth=2;
    h.Color = [0 0 0];
    set(gca,'ylim',yl);
end