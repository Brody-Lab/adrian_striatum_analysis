%% uses hacky script to add sessions to table that aren't sorted yet.

%% then overviews psychometrics by rat and makes a histogram of AP position by session.


T = load_sessions_table;

figure;
subplot(1,3,1);
h=histogram(T.craniotomy_AP,'NumBins',17,'Orientation','horizontal');

BinEdges = h.BinEdges;

rats = unique(T.rat);

for i=1:length(rats)
    aps{i} = T.craniotomy_AP(strcmp(T.rat,rats{i}));
    counts(:,i) = histcounts(aps{i},BinEdges);
end

colors = distinguishable_colors(length(rats));cla;hold off;
b = barh(0.5*(BinEdges(2:end) + BinEdges(1:end-1)),counts,'stacked');

for i=1:length(b)
    b(i).FaceColor = colors(i,:);
    b(i).DisplayName = rats{i};
end
legend(b,'Location','southeast');
xlabel('Number of Sessions')
ylabel({'mm anterior','to Bregma'},'Rotation',0,'HorizontalAlignment','right');


subplot(1,3,2);

for i=1:length(rats)
    count=0;
    dates = T.sess_date(strcmp(T.rat,rats{i}));
    for k=1:length(dates)
        sessid = bdata(['select sessid from sessions where ratname="' rats{i} '" and sessiondate="' datestr(dates(k),'YYYY-mm-DD') '" and hostname regexp "Rig20.*"']);        
        if length(sessid)==1
            count=count+1;
            sessids{i}(count)=sessid;
        else
            error('');
        end
    end
    protocol_data(i) = getProtocolData(rats{i},sessids{i},'PBups.*','getPeh',true);                
end

data=pbups_psych('','','protocolData',protocol_data,'plot',true,'xval','bupDiff','nPsychBins',10,'plotSpecial',false,'varDotSize',false,'subplot',false);
fithandles = cellfun(@(x)x.fitHandle,data);
%legend(fithandles,rats);
datahandles = cellfun(@(x)x.dataHandle,data);
for i=1:numel(datahandles)
    datahandles(i).Visible="off";
    fithandles(i).LineWidth=1.5;
end
set(gca,'xlim',[-40 40]);
title('');
pos=[0.06 0.57 0.14 0.26];
set(gcf,'units','normalized','position',pos,'color',[1 1 1]);



   

% clear aps counts b
% T = load_cells_table;
% T=T(T.is_in_dorsal_striatum,:);
% recording_names = T.recording_name;
% tmp=char(recording_names);
% rat=string(tmp(:,1:4));
% 
% figure;
% h=histogram(T.AP,'NumBins',12,'Orientation','horizontal');
% 
% BinEdges = h.BinEdges;
% 
% rats = unique(rat);
% 
% for i=1:length(rats)
%     aps{i} = T.AP(strcmp(rat,rats{i}));
%     counts(:,i) = histcounts(aps{i},BinEdges);
% end
% 
% colors = distinguishable_colors(length(rats));cla;hold off;
% b = barh(0.5*(BinEdges(2:end) + BinEdges(1:end-1)),counts,'stacked');
% 
% for i=1:length(b)
%     b(i).FaceColor = colors(i,:);
%     b(i).DisplayName = rats{i};
% end
% legend(b,'Location','southeast');
% xlabel('Number of Cells')
% ylabel({'mm anterior','to Bregma'},'Rotation',0,'HorizontalAlignment','right');
% set(gca,'xscale','linear');
% 
% clear aps counts b


clear aps counts b
T = load_cells_table;
T=T(T.is_in_dorsal_striatum,:);

subplot(1,3,3);
h=histogram(T.AP,'NumBins',25,'FaceColor','k');

ylabel('Number of Cells')
xlabel({'mm anterior','to Bregma'});
set(gca,'yscale','log','ytick',[10 100 1000],'yticklabel',[10 100 1000],'XDir','reverse');
jf=java.text.DecimalFormat; % comma for thousands, three decimal places
numOut= char(jf.format(height(T)));
text(0,4000,sprintf('n = %s neurons',numOut),'FontSize',15);

