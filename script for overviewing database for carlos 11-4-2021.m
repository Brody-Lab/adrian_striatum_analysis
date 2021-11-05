%% uses hacky script to add sessions to table that aren't sorted yet.

%% then overviews psychometrics by rat and makes a histogram of AP position by session.


T = add_missing_sessions_to_sessions_table();

figure;
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


figure;

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

pbups_psych('','','protocolData',protocol_data,'plot',true,'xval','bupDiff','nPsychBins',10,'plotSpecial',false,'varDotSize',false);



   