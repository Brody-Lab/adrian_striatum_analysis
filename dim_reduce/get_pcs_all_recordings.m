% get pcs for striatal populations in all recordings

min_dims = 10;
nsubsamples=100;

params = {'ref_event','cpoke_in',...
          'resolution',0.01,...
          'time_window_s',[0 3]...
%           'thresh_stability',[],...
%           'thresh_presence',[],...
%           'thresh_meanrate',[],...
%           'norm_factor',[],...
         };

paths = get_data_paths();
rng(1);
stats(numel(paths)).ncells=[];
parfor i=1:numel(paths)
    fprintf('\nLoading file %g of %g: %s \n',i,numel(paths),paths(i).cells_file);
    Cells = load(paths(i).cells_file);
    units = select_units(Cells,params{:},'exclude_cells',~Cells.is_in_dorsal_striatum);
    stats(i).ncells=sum(units);
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

for i=1:numel(paths)
    if stats(i).ncells>=min_dims
        count=count+1;
        idx(count)=i;
        var_explained(:,count) = mean([stats(i).pca_output.explained],2);
    end
end

cum_var_explained = cumsum(var_explained);


 mean(cum_var_explained(2,ap_group(idx)==1))
  mean(cum_var_explained(2,ap_group(idx)==2))
   mean(cum_var_explained(2,ap_group(idx)==3))
    mean(cum_var_explained(2,ap_group(idx)==4))
    
dims=[stats.ncells];    
figure;
scatter(dims(ap_group==4 & dims>=min_dims),cum_var_explained(2,ap_group(idx)==4)); hold on
scatter(dims(ap_group==3 & dims>=min_dims),cum_var_explained(2,ap_group(idx)==3))
scatter(dims(ap_group==2 & dims>=min_dims),cum_var_explained(2,ap_group(idx)==2))
scatter(dims(ap_group==1 & dims>=min_dims),cum_var_explained(2,ap_group(idx)==1))
xlabel('Number of Neurons')
ylabel({'% Variance Explained','by Top 5 PCs'})
