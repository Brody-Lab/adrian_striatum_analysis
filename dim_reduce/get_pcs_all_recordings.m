get pcs for striatal populations in all recordings

P=get_parameters;
min_dims = 20;
nsubsamples=100;

params = {'ref_event','cpoke_in',...
          'resolution',0.05,...
          'time_window_s',[0.5 1.5]...
           'thresh_stability',20,...
           'thresh_presence',0.01,...
           'thresh_meanrate',0,...
           'norm_factor',5,...
         };

paths = get_data_paths();
%rng(1);
parfor i=1:numel(paths)
    fprintf('\nLoading file %g of %g: %s \n',i,numel(paths),paths(i).cells_file);
    Cells = load(paths(i).cells_file);
    units = select_units(Cells,params{:},'exclude_cells',~Cells.is_in_dorsal_striatum);
    ncells(i)=sum(units);
    if sum(units)>=min_dims
        exclude_trials = validate_trials(Cells.Trials,'mode','agb_glm','quiet',true);
        if sum(units)==min_dims     
            stats(i) = get_pcs(Cells,'trial_idx',~exclude_trials,params{:},'units',units,'npcs',min_dims,'nsubsamples',1);            
        else
            tic;stats(i) = get_pcs(Cells,'trial_idx',~exclude_trials,params{:},'units',units,'npcs',min_dims,'nsubsamples',nsubsamples,'precision','single','gpu',false);toc
        end
        stats(i).pca_output = rmfield(stats(i).pca_output,'score'); % this takes up >90% of the memory and is not used at all in this script
    else
       fprintf('Skipping session because the number of selected units (%d) is less than the minimum (%d).\n',sum(units),min_dims); 
    end
    ap_group(i)=Cells.ap_group;
end

count=0;
clear idx
clear var_explained

for i=1:numel(paths)
    if stats(i).ncells>=min_dims
        count=count+1;
        idx(count)=i;
        var_explained(:,count) = mean([stats(i).pca_output.explained],2);
    end
end

cum_var_explained = cumsum(var_explained);


 mean(cum_var_explained(5,ap_group(idx)==1))
  mean(cum_var_explained(5,ap_group(idx)==2))
   mean(cum_var_explained(5,ap_group(idx)==3))
    mean(cum_var_explained(5,ap_group(idx)==4))
    
dims=ncells;
figure;
scatter(dims(ap_group==4 & dims>=min_dims),cum_var_explained(5,ap_group(idx)==4)); hold on
scatter(dims(ap_group==3 & dims>=min_dims),cum_var_explained(5,ap_group(idx)==3))
scatter(dims(ap_group==2 & dims>=min_dims),cum_var_explained(5,ap_group(idx)==2))
scatter(dims(ap_group==1 & dims>=min_dims),cum_var_explained(5,ap_group(idx)==1))
xlabel('Number of Neurons')
ylabel({'% Variance Explained','by Top 5 PCs'})



figure('color','w');
subplot(1,2,1);
for k=1:2
    for i=1:4
        boots = bootstrp(1000,@nanmean,cat(1,zeros(1,sum(ap_group(idx)==i)),cum_var_explained(:,ap_group(idx)==i))'./100);
        if k==1
            h(i)=shadedErrorBar(0:min_dims,boots,{@mean,@std},{'color',P.ap_group_colors(i,:)});hold on;
            h(i).patch.FaceAlpha=0.8;        
            h(i).mainLine.LineWidth=2;
        else
            plot(0:min_dims,mean(boots),'color',P.ap_group_colors(i,:),'LineWidth',3);hold on;
        end
    end
end
set(gca,P.axes_properties{:});
ylabel('Cumulative Fraction Variance Explained');
xlabel('Number of PCs');
set(gca,'xlim',[0 15]);
legend([h.mainLine],P.ap_group_labels);

subplot(1,2,2);

for i=1:size(cum_var_explained,2)
    dimensionality(i)=find(cum_var_explained(:,i)>=80,1);
end

for i=1:4
    boots=bootstrp(1000,@nanmean,dimensionality(ap_group(idx)==i));hold on
    errorbar(i,mean(boots),std(boots),'color',P.ap_group_colors(i,:),'CapSize',0,'Marker','.','MarkerSize',20,'LineWidth',2);
end
set(gca,'xtick',1:4,'xticklabel',P.ap_group_labels_xtick,P.axes_properties{:},'ytick',[5 6 7 8 9],'xlim',[0.5 4.5]);
ylabel('PCs needed to explain >80% variance');
