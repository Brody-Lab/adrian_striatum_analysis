% plotting gpfa latents scratch


% means
boots=bootstrp(1000,@mean,l);
boots_old=bootstrp(1000,@mean,l_old);
boots_std=bootstrp(1000,@std,l);
boots_old_std=bootstrp(1000,@std,l_old);
figure;
subplot(2,2,1);
b=boxplot(boots_old,'PlotStyle','compact','Notch','off','OutlierSize',1e-5,'MedianStyle','line');
set(gca,'xtick',1:5);
subplot(2,2,2);
b=boxplot(boots,'PlotStyle','compact','Notch','off','OutlierSize',1e-5,'MedianStyle','line');
set(gca,'xtick',1:5);
subplot(2,2,3);
b=boxplot(boots_old_std,'PlotStyle','compact','Notch','off','OutlierSize',1e-5,'MedianStyle','line');
set(gca,'xtick',1:5,'yscale','log','ylim',[1e-3 1e1]);
subplot(2,2,4);
b=boxplot(boots_std,'PlotStyle','compact','Notch','off','OutlierSize',1e-5,'MedianStyle','line');
set(gca,'xtick',1:5,'yscale','log','ylim',[1e-3 1e1]);



%% plotting click aligned latents (new and old set of latents)
figure;
for i=1:5
    subplot(1,2,1);
    hold on;plot(kPETH.timeS.clicks_on,gpfa_old.psth(i).clicks_on);
    title('Old (slower) latents')    
    xlabel('Time (s) after clicks on');
    set(gca,'ylim',[-3 3]);    
    subplot(1,2,2);
    hold on;plot(kPETH.timeS.clicks_on,gpfa.psth(i).clicks_on);    
    title('New (faster) latents')
    xlabel('Time (s) after clicks on');    
    set(gca,'ylim',[-3 3]);
end