
%% clustering
% run "glmfit_2022_03_03_10_23_17";
bs=[fits_table.stats.beta];
bs=bs(2:83,:);
ap_group = equalize_ap_groups(fits_table.ap_group);
bs=bs(:,~isnan(ap_group));
ap_group = ap_group(~isnan(ap_group));


responsiveness = max((bs.^2));
cutoff = prctile(responsiveness,40);

bs= bs(:,responsiveness>cutoff);
ap_group = ap_group(responsiveness>cutoff);


n=14;
% clustering approaches
T2 = clusterdata(bs','maxclust',n,'linkage',"ward");
T2 = cluster(linkage(bs',"ward"),'cutoff',1.5,'Criterion','Distance');
T = kmeans(bs',n);
%T = dbscan(bs',0.6,100);

% histograms of group membership
figure;for i=1:14;subplot(2,7,i);histogram(ap_group(T==i));set(gca,'xlim',[0.5 4.5]);end;matchylim(gcf);
figure;for i=1:4;subplot(1,4,i);histogram(T(ap_group==i));end



%% tsne
% euclidean and cosine both work well.
Y=tsne(bs','NumDimensions',2,'Distance','euclidean','Perplexity',30);

figure;
for i=1:max(T)
colors=jet(max(T));
hold on;scatter(Y(T==i,1),Y(T==i,2),25,colors(i,:),'Marker','.');
end


figure;
for i=1:4
colors=jet(4);
hold on;scatter(Y(ap_group==i,1),Y(ap_group==i,2),45,colors(i,:),'Marker','.');
end


hold on;
for i=1:4
colors=jet(4);
hold on;scatter(mean(Y(ap_group==i,1)),mean(Y(ap_group==i,2)),50,colors(i,:),'Marker','x');
end