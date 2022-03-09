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
    p.addParameter('reorder_cells',false,@(x)validateattributes(x,{'logical'},{'scalar'}));  
    p.addParameter('show_cell_no',false,@(x)validateattributes(x,{'logical'},{'scalar'}));            
    p.addParameter('smoothing_window_size',0,@(x)validateattributes(x,{'numeric'},{'nonnegative','nonempty','scalar'})); 
    p.addParameter('var_cutoff',1e-3,@(x)validateattributes(x,{'numeric'},{'positive','nonempty','scalar'}));        
    p.parse(varargin{:})
    params=p.Results;
    
    %% get weight matrix
    if ~params.use_combined_weights
        bs = [fits_table.stats.beta]';
        if params.biasCol
            bs=bs(:,2:end); % remove bias column
        end
        [column_idx,edim] = select_beta_columns(dspec,params.covariates);
        bs=bs(:,column_idx);
    else
        % get combined weight matrix
        stats=[fits_table.stats];
        for i=numel(params.covariates):-1:1
            [ws,tr{i}] = get_combined_weights_downsample(stats,params.covariates{i});
            bs{i} = ws.data;
            good_idx = var(bs{i})>params.var_cutoff;
            bs{i} = bs{i}(:,good_idx);
            tr{i} = tr{i}(:,good_idx);
            edim(i) = size(bs{i},2);
        end
        bs = cat(2,bs{:});
    end
    edim=cumsum(edim);    
    
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
    old_bs=bs;
    cv = cov(bs);
    mn = mean(bs);
    %bs=mvnrnd(mn(:),cv,size(bs,1));
    %bs=rand(size(bs));
    count=0;
    val=[3 4 5 6 8 10 12 14 16 18 20 30 50 100];
    for k=val
        count=count+1;
    for i=1:50
    n=numel(params.group_by);
    idx(:,i)= randsample(n,round(n*0.95),false);
    bs=old_bs(idx(:,i),:);
    %params.group_by = params.group_by(idx);
    ngroups= numel(unique(params.group_by));
    ncells = numel(params.group_by);
    
    %% determine tree structure
    tree = linkage(bs,params.linkage_method,params.metric);
    
    %% determine distance cutoff to produce desired number of clusters   
    clusterfun = @(cutoff)cluster(tree,'Cutoff',cutoff,'Criterion','distance'); 
    [cutoff,T(:,i)] = find_cutoff(clusterfun,k,5,10);
    end
    ar2(count) = get_RI(T,idx);
    end
    %T = spectralcluster(bs,params.max_clust);  
    %T=kmeans(bs,params.max_clust);
    %[~,outperm] = sort(T);
    %bs=bs(outperm,:);
    %params.group_by = params.group_by(outperm);      
    %end
    %% cluster and plot the color-thresholded dendrogram
    figure('Units','normalized','Position',[0.05 0.05 0.74 0.87],'color','w');
    dendrogram_ax = axes('Position',[0.03 0.1 0.07 0.8]);
    [~,~,outperm] = dendrogram(tree,Inf,'orientation','left','ColorThreshold',cutoff,'CheckCrossing',true);
    if params.reorder_cells
        [~,idx] = sort(params.group_by);
        bs = bs(idx,:);
    else
        bs=bs(outperm,:);
    end
    params.group_by = params.group_by(outperm);
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);axis off
    
    %% compute membership fractions
    clusters_in_order = unique(T(outperm),'stable');
    for i=1:params.max_clust
        these_groups = params.group_by(T(outperm)==clusters_in_order(i));
        n_cells_per_group(i)  = numel(these_groups);
        fracs(i,:) = histcounts(these_groups,[0.5 1.5 2.5 3.5 4.5])./n_cells_per_group(i);
    end
    
    %% if desired, reorder cells by cluster membership fraction
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
        cla(dendrogram_ax); % clear dendrogram since the cells are not in that order anymore
    elseif params.reorder_cells
        cla(dendrogram_ax);
    end
        
    
    %% optional smooth across cells for visualization
    if params.smoothing_window_size
        bs = filterArray(bs,gausswin(params.smoothing_window_size));
    end
    %% draw the image
    click_ticks = [0 0.1 1];
    covariate_ticks = [-1 -0.5 0 0.5 1];
    bs=log10(exp(bs));
    clim = [min(bs(:)) max(bs(:))]/2;
    clim=(clim-mean(clim));
    time_range = cellfun(@range,tr);
    gap=0.01;
    if params.reorder_cells
        linepos = find(diff(sort(params.group_by)));
    else
        linepos = cumsum(n_cells_per_group);            
    end    
    for i=1:numel(params.covariates)  
        if tr{i}(1)>0
            tr{i} = tr{i}-tr{i}(1);
        end
        if i==1
            these_cols = 1:edim(1);
            pos = [0.12 0.1 4e-2*time_range(i) 0.8];     % change this 4e-2 scaling value if you want to stretch the images more or less        
        else
            these_cols = (edim(i-1)+1) : edim(i);
            pos = [pos(1)+gap+pos(3) 0.1 4e-2*time_range(i) 0.8]; % and here. this is a bit hacky and i should automatically pick this value to cover the desired range of the images in the plot          
        end
        img_axes(i) = axes('Position',pos);                  
        if contains(params.covariates{i},"click")        
            imagesc(1:numel(tr{i}),1:ncells,bs(:,these_cols));set(gca,'clim',clim);colormap(redblue);
        else
            imagesc(tr{i},1:ncells,bs(:,these_cols));set(gca,'clim',clim);colormap(redblue);            
        end
        xl=get(gca,'xlim'); 
        yl=get(gca,'ylim'); 
        ax=gca;        
        if i==1
            xlabel('Time (s)');   
            if (params.group_reorder || params.reorder_cells) && params.show_cell_no
                set(gca,'ytick',[0:100:ncells]);h=ylabel('Cells');  
                yruler = ax.YRuler;
                yruler.Axle.Visible = 'off';            
            else
                ax.YAxis.Visible='off';                  
            end
        else
            ax.YAxis.Visible='off';  
        end
        if contains(params.covariates{i},"click")
            for k = 1:numel(click_ticks)
                [~,zr(k)] = min(abs(tr{i}-click_ticks(k)));
            end
            set(gca,'xtick',zr,'xticklabel',click_ticks,P.axes_properties{:}); box off;                 
            line(ones(1,2)*zr(1),yl,'LineWidth',1,'LineStyle','--','color',ones(1,3)/2);
        else
            set(gca,'xtick',covariate_ticks,P.axes_properties{:});box off;   
            line([0 0],yl,'LineWidth',1,'LineStyle','--','color',ones(1,3)/2);            
        end
        if params.reorder_cells
            for k=1:(ngroups-1)
                line(xl,ones(1,2)*linepos(k),'Color','k','LineWidth',1 );
            end            
        else          
            for k=1:params.max_clust
                line(xl,ones(1,2)*linepos(k),'Color','k','LineWidth',1 );
            end 
        end
        ax.FontSize = 12;           
        title(strsplit(P.covariate_names{i},' '));
        if i == numel(params.covariates)
            originalSize = get(gca, 'Position');            
            h=colorbar('Location','manual');
            h.Label.String = 'Gain';
            h.Label.Position = [0.5 2 0];
            h.Label.Rotation=0;
            h.FontSize=12;
            h.Box='off';
            h.LineWidth=0.75;
            h.TickLength=0.05;
            h.Position = [pos(1)+gap+pos(3) 0.8 0.01 0.1];
            h.Ticks = log10([0.01:0.01:0.09 0.1:0.1:1 2 3 4 5 6 7 8 9 10 20 30 40]);
            idx = ~ismember(h.Ticks,[-2 -1 0 1]);
            TickLabels = arrayfun(@num2str,[0.01:0.01:0.09 0.1:0.1:1 2 3 4 5 6 7 8 9 10 20 30 40],'uni',0);
            TickLabels(idx) = deal({''});
            h.TickLabels = TickLabels;
            set(gca,'Position',originalSize);
        end
    end


    if ~params.reorder_cells
        axes('Position',[0.85 0.1 0.1 0.8]) ;    
        h=barh(fracs,'stacked');axis off;xl=get(gca,'xlim');set(gca,'ydir','reverse','ylim',[0.5 params.max_clust+0.5],'xlim',[eps xl(2)]);
        for i=1:4;h(i).FaceColor = P.ap_group_colors(i,:);end    
        l=legend(h,P.ap_group_labels);
        l.FontSize=12;
        l.Position = [0.85,0.025,0.1,0.06];
    else
        groups = sort(params.group_by);        
        for i=1:ngroups
            height = mean(find(groups == i));
            axes(img_axes(1));
            text(-2,height,P.ap_group_labels_xtick{i},'color',P.ap_group_colors(i,:),'FontSize',16,'HorizontalAlignment','center');
        end       
    end
   

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