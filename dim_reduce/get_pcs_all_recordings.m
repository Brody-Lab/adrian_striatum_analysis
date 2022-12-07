% get pcs for striatal populations in all recordings


params = {'ref_event','cpoke_in',...
          'resolution',0.01,...
          'time_window_s',[0 3]...
%           'thresh_stability',[],...
%           'thresh_presence',[],...
%           'thresh_meanrate',[],...
%           'norm_factor',[],...
         };

paths = get_data_paths();

parfor i=1:numel(paths)
    fprintf('\nLoading file %g of %g: %s \n',i,numel(paths),paths(i).cells_file);
    Cells = load(paths(i).cells_file);
    units = select_units(Cells,params{:},'exclude_cells',~Cells.is_in_dorsal_striatum);
    if any(units)
        exclude_trials = validate_trials(Cells.Trials,'mode','agb_glm','quiet',true);
        stats(i) = get_pcs(Cells,'trial_idx',~exclude_trials,params{:},'units',units);
    else
       fprintf('Skipping session because there are no selected units.\n'); 
       stats(i).dim=NaN;
    end
    ap_group(i)=Cells.ap_group;
end

dims=[stats.dim];
min_dims = 15;
max_dims = max(dims);
var_explained = NaN(max_dims,sum(dims>=min_dims));
count=0;
clear idx

for i=1:numel(paths)
    if dims(i)>=min_dims
        count=count+1;
        idx(count)=i;
        var_explained(1:dims(i),count) = stats(i).explained;
    end
end

cum_var_explained = cumsum(var_explained);

norm_cum_var_explained=cum_var_explained ./ cum_var_explained(min_dims,:);

 mean(norm_cum_var_explained(5,ap_group(idx)==1))
  mean(norm_cum_var_explained(5,ap_group(idx)==2))
   mean(norm_cum_var_explained(5,ap_group(idx)==3))
    mean(norm_cum_var_explained(5,ap_group(idx)==4))
    
    
figure;
scatter(dims(ap_group==4 & dims>=min_dims),cum_var_explained(5,ap_group(idx)==4)); hold on
scatter(dims(ap_group==3 & dims>=min_dims),cum_var_explained(5,ap_group(idx)==3))
scatter(dims(ap_group==2 & dims>=min_dims),cum_var_explained(5,ap_group(idx)==2))
scatter(dims(ap_group==1 & dims>=min_dims),cum_var_explained(5,ap_group(idx)==1))
xlabel('Number of Neurons')
ylabel({'% Variance Explained','by Top 5 PCs'})
