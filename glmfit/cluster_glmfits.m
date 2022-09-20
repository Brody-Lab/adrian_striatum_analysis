function cluster_glmfits(fits_table,dspec,varargin)

    % takes a set of glm fits and displays them as an image where rows are
    % cells and clumns are covariates. optionally perform hierarchical
    % clustering to determine cell ordering. 

    %% parse and validate inputs
    P=get_parameters();
    p=inputParser;
    default_covariates = P.covariate_order(ismember(P.covariate_order,{dspec.covar.label}));
    p.addParameter('covariates',default_covariates,@(x)validateattributes(x,{'cell'},{'nonempty'}));
    p.addParameter('covariate_names',P.covariate_names(ismember(P.covariate_order,{dspec.covar.label})),@(x)validateattributes(x,{'cell'},{'nonempty'}));   
    p.addParameter('covariates_to_plot',default_covariates,@(x)validateattributes(x,{'cell'},{'nonempty'}));
    p.addParameter('covariates_to_cluster',default_covariates,@(x)validateattributes(x,{'cell'},{'nonempty'}));    
    p.addParameter('max_clust',[],@(x)validateattributes(x,{'numeric'},{'positive','nonempty','scalar'}));
    p.addParameter('linkage_method','ward',@(x)validateattributes(x,{'char','string'},{'nonempty'}));
    p.addParameter('metric','euclidean',@(x)validateattributes(x,{'char','string'},{'nonempty'}));
    p.addParameter('equalize_ap_groups',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('responsive_fun',@(x)max(x.^2),@(x)validateattributes(x,{'function_handle'},{}));
    p.addParameter('responsive_cutoff_prctile',0,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',100}));
    p.addParameter('biasCol',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('use_combined_weights',true,@(x)validateattributes(x,{'logical'},{'scalar'})); % if false, things will break. I've found using the combined weights so much better that code assumes it in some places.
    p.addParameter('group_by',ones(size(fits_table,1),1)); % should include ascending natural numbers starting at 1 (or NaNs to be excluded).
    p.addParameter('group_reorder',true,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('reorder_cells',false,@(x)validateattributes(x,{'logical'},{'scalar'}));  
    p.addParameter('show_cell_no',false,@(x)validateattributes(x,{'logical'},{'scalar'}));      
    p.addParameter('find_optimal_n_clust',false,@(x)validateattributes(x,{'logical'},{'scalar'}));     
    p.addParameter('variable_bar_size',false,@(x)validateattributes(x,{'logical'},{'scalar'}));                  
    p.addParameter('smoothing_window_size',0,@(x)validateattributes(x,{'numeric'},{'nonnegative','nonempty','scalar'})); 
    p.addParameter('var_cutoff',1e-3,@(x)validateattributes(x,{'numeric'},{'positive','nonempty','scalar'}));       
    p.addParameter('clim_multiplier',5);
    p.addParameter('group_colors',P.ap_group_colors);
    p.addParameter('group_labels',P.ap_group_labels);
    p.addParameter('click_ticks',[0 0.1 1]);
    p.addParameter('covariate_ticks',-5:0.5:5);    % ticks outside of the plotted range are ignored
    p.addParameter('covariate_gap',0.01); % horizontal gap between covariates
    p.addParameter('link','log');
    p.addParameter('colormap',redblue);
    p.addParameter('fig_pos',[0.05 0.05 0.6 0.8]);
    p.parse(varargin{:})
    params=p.Results;
        
    %% get weight matrix
    [bs,edim,tr] = get_data_matrix(fits_table.stats,dspec,params);
    bs(isinf(bs))=NaN;
    
    %% remove bad elements (where group variable is NaN or unresponsive)
    [bs,params.group_by] = remove_bad_rows(bs,params);
    params.ngroups= numel(unique(params.group_by));
    params.ncells = numel(params.group_by);    
    
    %% get columns to use for clustering (if you want clustering to be performed on only certain covariates)
    clust_idx = find(ismember(params.covariates,params.covariates_to_cluster));
    e  = [0 edim];
    count=0;
    for i=clust_idx(:)'
        count=count+1;
        clust_cols{count} = [(e(i)+1):e(i+1)]';
    end
    match_covariate_cols=true;
    if match_covariate_cols
        count=0;
        min_n=min(cellfun(@numel,clust_cols));
        for i=clust_idx(:)'
            count=count+1;            
            clust_cols{count} = randsample(clust_cols{count},min_n,false);
        end
    end
    clust_cols=cat(1,clust_cols{:});
    
    %% optionally find optimal number of clusters
    if params.find_optimal_n_clust   
        [params.max_clust,params.ARI,params.ARI_rand,params.vals] = find_optimal_n_clust(bs(:,clust_cols));
    end
    
    %% determine tree structure
    tree = linkage(bs(:,clust_cols),params.linkage_method,params.metric);
    
    %% determine distance cutoff to produce desired number of clusters   
    clusterfun = @(cutoff)cluster(tree,'Cutoff',cutoff,'Criterion','distance'); 
    [cutoff,T] = find_cutoff(clusterfun,params.max_clust,5,10);

    %% initialize figure
    figure('Units','normalized','Position',params.fig_pos,'color',[1 1 1]);    
    
    %% compute dendrogram order
    [~,~,outperm] = dendrogram(tree,0,'orientation','left','ColorThreshold',cutoff,'CheckCrossing',false);
    if params.reorder_cells || params.group_reorder
        clf;
    else
        params.dendrogram_pos = [0.03 0.1 0.07 0.8];        
        set(gca,'position',params.dendrogram_pos,'ydir','reverse','xtick',[],'ytick',[]);axis off 
    end
    
    %% rearrange according to grouping variable if desired
    if params.reorder_cells
        [~,outperm] = sort(params.group_by);
    end    
    
    %% resort by desired order
    bs=bs(outperm,:);
    params.group_by = params.group_by(outperm);        
    
    %% if desired, reorder cells by cluster membership fraction (group reorder)
    if params.group_reorder
        clusters_in_order = unique(T(outperm),'stable');
        for i=params.max_clust:-1:1
            these_groups = params.group_by(T(outperm)==clusters_in_order(i));
            n_cells_per_group(i)  = numel(these_groups);
            params.fracs(i,:) = histcounts(these_groups,(0:params.ngroups) + 0.5)./n_cells_per_group(i);
        end        
        avg_group = wmean(repmat(1:params.ngroups,params.max_clust,1),params.fracs,2);   
        [~,idx] = sort(avg_group);
        bs_old = bs;
        group_idx=[];
        bs=[];
        for i=1:params.max_clust
            group_idx = [group_idx;params.group_by(T(outperm) == clusters_in_order(idx(i)))];
            bs = [bs;bs_old(T(outperm) == clusters_in_order(idx(i)),:)];
        end
        params.fracs = params.fracs(idx,:);       
        n_cells_per_group = n_cells_per_group(idx);
    end
        
    %% optional smooth across cells for visualization
    if params.smoothing_window_size
        bs = filterArray(bs,gausswin(params.smoothing_window_size));
    end
    
    %% only keep covariates being plotted
    [bs,tr,edim] = select_covariates_to_plot(bs,tr,edim,params) ;  
    
    %% setting up the plotting
    params.clim  = [-1 1] * params.clim_multiplier * std(bs(:),'omitnan'); % clim is set relative to the s.d. to be insensitive to outliers
    time_range = cellfun(@range,tr);
    if params.reorder_cells
        params.linepos = find(diff(sort(params.group_by)));
    else
        params.linepos = cumsum(n_cells_per_group);            
    end    
    params.width = 0.82 - 0.115; % this defines the edges of the main plotting area (covariates), with dendrogram optionally on the left and group fractions optionally on the right
    params.effective_width = params.width - (numel(params.covariates_to_plot)-1)*params.covariate_gap; % this defines the total time range, i.e. width minus sum of the gaps.
    w = params.effective_width*time_range./sum(time_range); % width of each covariate, scaled by its time range   
    
    %% loop over covariates and plot
    for i=1:numel(params.covariates_to_plot)  
        
        % get data matrix indices and axis position
        if i==1
            these_cols = 1:edim(1);
            pos = [0.115 0.1 w(i) 0.8];                 
        else
            these_cols = (edim(i-1)+1) : edim(i);
            pos = [0.115 + params.covariate_gap*(i-1)+sum(w(1:i-1)) 0.1 w(i) 0.8];             
        end   
        
        % main plotting routine
        plot_this_covariate(params,tr{i},bs(:,these_cols),pos,params.covariates_to_plot{i},params.covariate_names{i});
        
        % add cell numbers, group labels and time label if first covariate
        if i==1
            if params.reorder_cells
                groups = sort(params.group_by);        
                for k=1:params.ngroups
                    height = mean(find(groups == k));
                    text(-2.4,height,params.group_labels{k},'color',params.group_colors(k,:),'FontSize',22,'HorizontalAlignment','center');
                end    
            end  
            g=xlabel('Time (s)');  
            g.FontSize=20;
            if (params.group_reorder || params.reorder_cells) && params.show_cell_no
                set(gca,'ytick',0:100:params.ncells);
                ylabel('Cells');  
                yruler = ax.YRuler;
                yruler.Axle.Visible = 'off';            
            else
                ax.YAxis.Visible='off';                  
            end
        end
    end
    
    %% make colorbar
    draw_colorbar(pos,params);
        
    %% set colormap
    colormap(params.colormap);

    %% draw bars indicating fraction of clusters from each group 
    draw_cluster_fractions(params);
 
end

function [bs,edim,tr] = get_data_matrix(stats,dspec,params)
    if ~params.use_combined_weights
        bs = [stats.beta]';
        if params.biasCol
            bs=bs(:,2:end); % remove bias column
        end
        [column_idx,edim] = select_beta_columns(dspec,params.covariates);
        bs=bs(:,column_idx);
        tr=[];
    else
        % get combined weight matrix
        for i=numel(params.covariates):-1:1
            [ws,tr{i}] = get_combined_weights_downsample(stats,params.covariates{i},10); % downsample by a factor of 10
            bs{i} = ws.data;
            good_idx = var(bs{i})>params.var_cutoff;
            bs{i} = bs{i}(:,good_idx);
            tr{i} = tr{i}(:,good_idx);
            if tr{i}(1)>0 % if the first xval is slightly above zero, just set it to zero so we get nicer ticks.
                tr{i} = tr{i}-tr{i}(1);
            end            
            edim(i) = size(bs{i},2);
        end
        bs = cat(2,bs{:});
    end
    edim=cumsum(edim);    
    if params.link=="log"
        bs=log10(exp(bs)); % convert to log base 10 for interpretability
    else
        error(''); % not implemented yet because not clear what you want to do here unless given an example
    end
end

function [bs,group_by] = remove_bad_rows(bs,params)
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
    
    group_by = params.group_by;
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
    max_count=100;
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
            %fprintf('Found cutoff of %g for %g clusters in %g iterations.\n',cutoff,n_clust,count);
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
        if count>max_count
            cutoff=guess;
            fprintf('Could not find cutoff for %g clusters in %g iterations. Made %g instead.\n',max_clust,count,n_clust);            
            return
        end
    end
end

function AR2 = get_RI(T,idx)

count=0;
for i=1:size(T,2)
    for k=1:size(T,2)
        if k<=i
            continue
        else
            count=count+1;
            [~,ia,ib] = intersect(idx(:,i),idx(:,k));
            [AR2(count),RI2(count),MI2(count),HI2(count)] = RandIndex(T(ia,i),T(ib,k));
        end
    end
end
AR2 = mean(AR2);
end

function [nclust,ARI,ARI_rand,nclusts] = find_optimal_n_clust(X)
    P=get_parameters;
    mn = mean(X);
    n = size(X,1);
    X_rand=mvnrnd(mn(:),cov(X),n); % generate a simulated dataset of the same size with the same first and second moments (i.e. mvnormal distribution)
    nrep=70; % resample 95% of the real and simulated data this many times for each cluster size
    count=0;
    nclusts=[3:30 35 40 45 50 60 70 80 90 100];
    for k=nclusts
        count=count+1;
        for i=nrep:-1:1
            if count==1
                idx(:,i)= randsample(n,round(n*0.95),false);  
            end
            T(:,i) = clusterdata(X(idx(:,i),:),'maxclust',k,'linkage',"ward",'criterion','distance','savememory','off');   
            T_rand(:,i) = clusterdata(X_rand(idx(:,i),:),'maxclust',k,'linkage',"ward",'criterion','distance','savememory','off');                        
        end
        ARI(count) = get_RI(T,idx);
        ARI_rand(count) = get_RI(T_rand,idx);
    end  
    [~,idx] = max(ARI-ARI_rand);
    nclust = nclusts(idx);
    
    % plot
    figure('color',[1 1 1],'units','normalized','position',[0.17 0.64 0.19 0.26]);
    plot(nclusts,ARI-ARI_rand,'k','LineWidth',1.5);
    yl=get(gca,'ylim');
    set(gca,'box','off',P.axes_properties{:},'FontSize',16,...
        'XTick',20:20:100,'xlim',[0 100],'ylim',[0 yl(2)]);
    ylabel({'Excess Cluster Stability','(ARI)'});
    xlabel('Number of Clusters');
    [~,max_idx] = max(ARI-ARI_rand);
    hold on;scatter(nclust,ARI(max_idx)-ARI_rand(max_idx),500,[0 1 0],'.');
    text(nclust+2,ARI(max_idx)-ARI_rand(max_idx)+0.01,sprintf('n=%g',nclust),'FontSize',16,'color','g');
    line(ones(1,2)*nclust,get(gca,'ylim'),'color','g','LineWidth',1.5,'LineStyle','--');    
end

function [bs,tr,edim] = select_covariates_to_plot(bs,tr,edim,params)
    bad_cols=[];
    for i=numel(params.covariates):-1:1
        if ~ismember(params.covariates{i},params.covariates_to_plot)
            if i==1
                these_cols = 1:edim(1);
            else
                these_cols = (edim(i-1)+1) : edim(i);
                
            end 
            bad_idx(i)=true;
            bad_cols = union(bad_cols,these_cols);
        else
            bad_idx(i)=false;
        end
    end
    edim = [edim(1) diff(edim)];
    edim = edim(~bad_idx);
    edim = cumsum(edim);
    bs(:,bad_cols)=[];
    tr=tr(~bad_idx);
end

function plot_this_covariate(params,tr,bs,pos,covariate_to_plot,covariate_name)
    P = get_parameters();
    axes('Position',pos,'ydir','reverse');     
    ax=gca;
    if contains(covariate_to_plot,"click")        
        imagesc(1:numel(tr),1:params.ncells,bs);
        for k = numel(params.click_ticks):-1:1
            [~,zr(k)] = min(abs(tr-params.click_ticks(k)));
        end        
        set(ax,'xlim',[1 numel(tr)],'xtick',zr,'xticklabel',params.click_ticks);             
        line(ones(1,2)*zr(1),ax.YLim,'LineWidth',1.5,'LineStyle','--','color',ones(1,3)/2);        
    else
        imagesc(tr,1:params.ncells,bs);
        set(ax,'xlim',[min(tr) max(tr)],'xtick',params.covariate_ticks);
        line([0 0],ax.YLim,'LineWidth',1.5,'LineStyle','--','color',ones(1,3)/2);           
    end
    for k=1:numel(params.linepos)
        line(ax.XLim,ones(1,2)*params.linepos(k),'Color','k','LineWidth',1.5 );
    end   
    set(ax,'FontSize',16,'clim',params.clim,'ylim',[0.5 params.ncells],P.axes_properties{:});box off;
    ax.YAxis.Visible='off';      
    title(strsplit(covariate_name,' '),'FontSize',20);
end

function draw_colorbar(pos,params)
    originalSize = get(gca, 'Position'); % remember current axis position so it can be reset in case colorbar squeezes it        
    h=colorbar('Location','manual');
    ticks = log10([0.01:0.01:0.09 0.1:0.1:1 2 3 4 5 6 7 8 9 10 20 30 40]);
    idx = ~ismember(ticks,[-2 -1 0 1]);
    TickLabels = arrayfun(@num2str,[0.01:0.01:0.09 0.1:0.1:1 2 3 4 5 6 7 8 9 10 20 30 40],'uni',0);
    TickLabels(idx) = deal({''});    
    set(h,'FontSize',14,'box','off','LineWidth',0.75,'TickLength',0.05,'Position',[pos(1)+params.covariate_gap+pos(3) 0.8 0.01 0.1],'Ticks',ticks,'TickLabels',TickLabels);
    set(h.Label,'String','Gain','Position',[0.5 2 0],'Rotation',0);    
    set(gca,'Position',originalSize);
end

function draw_cluster_fractions(params)
    if ~params.reorder_cells 
        if params.variable_bar_size
            for i=1:numel(params.linepos)
                if i==1
                    real_width = 0.8*params.linepos(i)./params.ncells;
                else
                    real_width = 0.8*diff(params.linepos([i-1 i]))./params.ncells;           
                end
                if real_width>0.02
                    width = real_width - 0.01;
                else
                    width = real_width;
                end
                axes('Position',[0.87  0.9-0.8*params.linepos(i)./params.ncells + (real_width-width)/2  0.1 width]);                                   
                h=barh(repmat(params.fracs(i,:),2,1),'stacked','BarWidth',1);
                axis off;
                xl=get(gca,'xlim');
                set(gca,'ydir','reverse','xlim',[eps xl(2)],'ylim',[0.5 1.5]);                    
                for k=1:params.ngroups
                    h(k).FaceColor = params.group_colors(k,:);
                    h(k).LineWidth=1;
                end                      
            end              
        else    
            axes('Position',[0.87 0.1 0.1 0.8]) ;    
            h=barh(params.fracs,'stacked');
            axis off;
            xl=get(gca,'xlim');
            set(gca,'ydir','reverse','ylim',[0.5 params.max_clust+0.5],'xlim',[eps xl(2)]);
            for i=1:4;h(i).FaceColor = params.group_colors(i,:);h(i).LineWidth=1;end    
        end   
        l=legend(h,params.group_labels);
        l.FontSize=14;
        l.Position = [0.87,0.025,0.1,0.06];           
    end
end