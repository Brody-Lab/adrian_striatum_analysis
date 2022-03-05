
%% clustering
% run "glmfit_2022_03_03_10_23_17";
bs=[fits_table.stats.beta];
bs=bs(2:83,:);
ap_group = equalize_ap_groups(fits_table.ap_group);
bs=bs(:,~isnan(ap_group));
n=10;
% clustering approaches
T = clusterdata(bs',"linkage","ward",'maxclust',n);
T = kmeans(bs',n);


% histograms of group membership
figure;for i=1:10;subplot(2,5,i);histogram(g(idx==i));set(gca,'xlim',[0.5 4.5]);end
figure;for i=1:4;subplot(1,4,i);histogram(T(g==i));end



%% tsne
% euclidean and cosine both work well.
Y=tsne(bs','NumDimensions',3,'NumPCAComponents',50,'Distance','cosine','Perplexity',18 );
colors=jet(4);figure;for i=1:4;hold on;scatter3(Y(g==i,1),Y(g==i,2),Y(g==i,3),25,colors(i,:),'Marker','.');end