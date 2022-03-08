function cluster_glmfits(fits_table,dspec,varargin)

    % run "glmfit_2022_03_03_10_23_17" used for testing this code


    %% parse and validate inputs
    P=get_parameters();
    p=inputParser;
    p.addParameter('covariates',P.covariate_order(ismember(P.covariate_order,{dspec.covar.label})),@(x)validateattributes(x,{'cell'},{'nonempty'}));
    p.addParameter('max_clust',14,@(x)validateattributes(x,{'numeric'},{'positive','nonempty','scalar'}));
    p.addParameter('linkage_method','ward',@(x)validateattributes(x,{'char','string'},{'nonempty'}));
    p.addParameter('metric','euclidean',@(x)validateattributes(x,{'char','string'},{'nonempty'}));
    p.addParameter('equalize_ap_groups',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('responsive_fun',@(x)max(x.^2),@(x)validateattributes(x,{'function_handle'},{}));
    p.addParameter('responsive_cutoff_prctile',0,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',100}));
    p.addParameter('biasCol',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('use_combined_weights',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('group_by',fits_table.ap_group);
    p.addParameter('group_reorder',true,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.parse(varargin{:});
    params=p.Results;
    
    %% get weight matrix
    if ~params.use_combined_weights
        bs = [fits_table.stats.beta]';
        if params.biasCol
            bs=bs(:,2:end); % remove bias column
        end
        [column_idx,edim] = select_beta_columns(dspec,params.covariates);
        bs=bs(:,column_idx);
        edim=cumsum(edim);        
    else
        % get combined weight matrix
        stats=[fits_table.stats];
        for i=1:numel(params.covariates)
            ws = get_combined_weights_downsample(stats,params.covariates{i});
            bs{i} = ws.data;
            edim(i) = size(bs{i},2);
        end
        bs = cat(2,bs{:});
        vars  =var(bs);
        good_idx = vars>1e-2;       
        bs=bs(:,good_idx);
        edim=cumsum(edim);
        edim=arrayfun(@(x)sum(good_idx(1:x)),edim);
    end
    
    %% remove elements where group variable is nan
    bs=bs(~isnan(params.group_by),:);
    params.group_by = params.group_by(~isnan(params.group_by));

    %% remove unresponsive cells (determined for now by amplitudes of fit coefficients)
    if params.responsive_cutoff_prctile>0
        responsiveness = params.responsive_fun(bs');
        cutoff = prctile(responsiveness,params.responsive_cutoff_prctile);        
        bs= bs(responsiveness>cutoff,:);
        params.group_by = params.group_by(responsiveness>cutoff);   
    end
    ngroups= numel(unique(params.group_by));
    ncells = numel(params.group_by);
    
    %% determine tree structure
    tree = linkage(bs,params.linkage_method,params.metric);
    
    %% determine distance cutoff to produce desired number of clusters
    
    clusterfun = @(cutoff)cluster(tree,'Cutoff',cutoff,'Criterion','distance'); 
    [cutoff,T] = find_cutoff(clusterfun,params.max_clust,5,10);

       
    %% cluster and plot the color-thresholded dendrogram
    figure;
    dendrogram_ax = axes('Position',[0.1 0.1 0.1 0.8]);
    [h,t,outperm] = dendrogram(tree,Inf,'orientation','left','ColorThreshold',cutoff,'CheckCrossing',true);
    bs=bs(outperm,:);
    params.group_by = params.group_by(outperm);
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);axis off
    axes('Position',[0.22 0.1 0.6 0.8]);  
    
    %% compute membership fractions
    clusters_in_order = unique(T(outperm),'stable');
    for i=1:params.max_clust
        these_groups = params.group_by(T(outperm)==clusters_in_order(i));
        n_cells_per_group(i)  = numel(these_groups);
        fracs(i,:) = histcounts(these_groups,[0.5 1.5 2.5 3.5 4.5])./n_cells_per_group(i);
    end
    
    %% if desired, reorder cells
    if params.group_reorder
        avg_group = wmean(repmat([1 2 3 4],params.max_clust,1),fracs,2);   
        [~,idx] = sort(avg_group);
        bs_old = bs;
        group_idx=[];
        bs=[];
        for i=1:params.max_clust
            group_idx = [group_idx;params.group_by(T(outperm) == clusters_in_order(idx(i)))];
            bs = [bs;bs_old(T(outperm) == clusters_in_order(idx(i)),:)];
        end
        fracs = fracs(idx,:);       
        n_cells_per_group = n_cells_per_group(idx);
        cla(dendrogram_ax);
    end
        
    
    %% slightly smooth across cells for visualization
    %filter = gausswin(2);
    %img_filt = filterArray(bs,filter);

    %

    %% draw the image
    imagesc(bs);a=get(gca,'clim');a=a-mean(a);set(gca,'clim',a/3);colormap(redblue);
    set(gca,'xtick',[],'ytick',[]);axis off
    xl=get(gca,'xlim'); 
    yl=get(gca,'ylim');
    linepos = cumsum(n_cells_per_group);
    for i=1:params.max_clust
        line(xl,ones(1,2)*linepos(i),'Color','k','LineWidth',1 );
    end
    for i=1:(numel(params.covariates)-1)
        if ismember(params.covariates{i},{'stereo_click','left_clicks','cpoke_out_left','spoke_left_miss','spoke_left_hit'})
            line(ones(1,2)*edim(i),yl,'Color',ones(1,3)/2,'LineWidth',0.5,'LineStyle','--');
        else
            line(ones(1,2)*edim(i),yl,'Color','k','LineWidth',1);            
        end
    end

    
    axes('Position',[0.84 0.1 0.1 0.8]) ;    
    h=barh(fracs,'stacked');axis off;xl=get(gca,'xlim');set(gca,'ydir','reverse','ylim',[0.5 params.max_clust+0.5],'xlim',[eps xl(2)]);
    for i=1:4;h(i).FaceColor = P.ap_group_colors(i,:);end    
   

end


function [column_idx,edim] = select_beta_columns(dspec,covariates)
    column_idx = buildGLM.getDesignMatrixColIndices(dspec, covariates);
    edim = cellfun(@numel,column_idx);
    column_idx = cat(1,column_idx{:});
end


function [cutoff,T] = find_cutoff(clusterfun,max_clust,guess,step)
    last_guess=guess;
    min_guess = eps;
    count=0;
    while true
        T = clusterfun(guess);
        n_clust = numel(unique(T));
        count=count+1;
        if n_clust>max_clust            
            if guess<last_guess
                step=step/2;
            end
            last_guess=guess;            
            guess=guess+step;
        elseif n_clust==max_clust
            cutoff=guess;
            fprintf('Found cutoff of %g for %g clusters in %g iterations.\n',cutoff,n_clust,count);
            return
        elseif n_clust<max_clust
            if guess>last_guess
                step=step/2;
            end
            last_guess=guess;            
            guess=guess-step;
        end
        if guess<min_guess
            guess=eps;
        end        
    end
end