
function plot_glmfits(fits_table,dspec,varargin)

    % takes a set of glm fits and displays them as an image where rows are
    % cells and clumns are covariates. optionally perform hierarchical
    % clustering to determine cell ordering. 

    %% parse and validate inputs
    P=get_parameters();
    p=inputParser;
    default_covariates = P.covariate_order(ismember(P.covariate_order,fieldnames(fits_table.stats(1).ws)));
    p.addParameter('covariates',default_covariates,@(x)validateattributes(x,{'cell'},{'nonempty'}));
    p.addParameter('covariate_names',P.covariate_names(ismember(P.covariate_order,{dspec.covar.label})),@(x)validateattributes(x,{'cell'},{'nonempty'}));   
    p.addParameter('covariates_to_plot',default_covariates,@(x)validateattributes(x,{'cell'},{'nonempty'}));
    p.addParameter('equalize_groups',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('responsive_fun',@(x)max(x.^2),@(x)validateattributes(x,{'function_handle'},{}));
    p.addParameter('responsive_cutoff_prctile',0,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',100}));
    p.addParameter('biasCol',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('group_by',ones(size(fits_table,1),1)); % should include ascending natural numbers starting at 1 (or NaNs to be excluded).
    p.addParameter('show_cell_no',false,@(x)validateattributes(x,{'logical'},{'scalar'}));      
    p.addParameter('smoothing_window_size',0,@(x)validateattributes(x,{'numeric'},{'nonnegative','nonempty','scalar'})); 
    p.addParameter('var_cutoff',1e-8,@(x)validateattributes(x,{'numeric'},{'positive','nonempty','scalar'}));       
    p.addParameter('clim_multiplier',5);
    p.addParameter('linkage_method','ward',@(x)validateattributes(x,{'char','string'},{'nonempty'}));
    p.addParameter('metric','euclidean',@(x)validateattributes(x,{'char','string'},{'nonempty'}));  
    p.addParameter('max_clust',5,@(x)validateattributes(x,{'numeric'},{'positive','nonempty','scalar'}));    
    p.addParameter('sig',[]);
    p.addParameter('group_colors',P.ap_group_colors);
    p.addParameter('group_labels',P.ap_group_labels_xtick);
    p.addParameter('click_ticks',[0 0.1 1]);
    p.addParameter('covariate_ticks',-5:0.5:5);    % ticks outside of the plotted range are ignored
    p.addParameter('covariate_gap',0.015); % horizontal gap between covariates
    p.addParameter('link','log');
    p.addParameter('colormap',redblue);
    p.addParameter('make_fig',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.parse(varargin{:})
    params=p.Results;
        
    %% get weight matrix
    [bs,edim,tr] = get_data_matrix(fits_table.stats,params);
    bs(isinf(bs))=NaN;
    
    %% remove bad elements (where group variable is NaN or unresponsive)
    [bs,params.group_by] = remove_bad_rows(bs,params);
    params.ngroups= numel(unique(params.group_by));
    params.ncells = numel(params.group_by);    
   
    %% make sure in ap order
    [~,outperm] = sort(params.group_by);    
    bs=bs(outperm,:);
    params.group_by = params.group_by(outperm);        
        
    %% optional smooth across cells for visualization
    if params.smoothing_window_size
        bs = filterArray(bs,gausswin(params.smoothing_window_size));
    end
    
    %% only keep covariates being plotted
    [bs,tr,edim] = select_covariates_to_plot(bs,tr,edim,params) ;  
    order=[];
    clear outperm
    for i=1:4
        idx=find(params.group_by==i);
        tree = linkage(bs(idx,:),params.linkage_method,params.metric);
        clusterfun = @(cutoff)cluster(tree,'Cutoff',cutoff,'Criterion','distance'); 
        [cutoff,T{i}] = find_cutoff(clusterfun,params.max_clust,5,10);
        [~,~,outperm{i}] = dendrogram(tree,0,'orientation','left','ColorThreshold',cutoff,'CheckCrossing',false);    cla;    
        b=[];
        for k=1:params.max_clust
           ggg(k)=mean(find(T{i}(outperm{i}) == k)); 
        end
        [~,ii] = sort(ggg) ;     
        for k=ii(:)'
           g=find(T{i} == k); 
           n=20;
           if numel(g)>=n
                g=randsample(g,n,false);
           else
               g=g;
           end
            o= arrayfun(@(x)find(outperm{i}==x),g);
            [~,th] = sort(o);
            g=g(th);           
           b=[b;g(:)];
        end
        nn(i)=numel(b);
        order = [order;idx(b)];
    end
    bs=bs(order,:); 
    params.group_by = params.group_by(order);
    
    
    if params.equalize_groups
        params.group_by = equalize_groups(params.group_by);
        bs=bs(~isnan(params.group_by),:);
        params.group_by = params.group_by(~isnan(params.group_by));
        for i=1:4
            nn(i)=sum(params.group_by==i);
        end
    end
    params.ncells =size(bs,1);

    
    
    
    clf;
    %% initialize figure
    if params.make_fig
        figure('Units','normalized','Position',params.fig_pos,'color',[1 1 1]);    
    end
   
    
    %% setting up the plotting
    params.clim  = [-1 1] * params.clim_multiplier * nanstd(bs(:)); % clim is set relative to the s.d. to be insensitive to outliers
    time_range = cellfun(@range,tr);
    params.linepos = cumsum(nn);
    params.width = 0.5; % this defines the edges of the main plotting area (covariates), with dendrogram optionally on the left and group fractions optionally on the right
    params.effective_width = params.width - (numel(params.covariates_to_plot)-1)*params.covariate_gap; % this defines the total time range, i.e. width minus sum of the gaps.
    w = params.effective_width*time_range./sum(time_range); % width of each covariate, scaled by its time range   
    
    %% loop over covariates and plot
    for i=1:numel(params.covariates_to_plot)  
        
        % get data matrix indices and axis position
        if i==1
            these_cols = 1:edim(1);
            pos = [0.1 0.1 w(i) 0.85];                 
        else
            these_cols = (edim(i-1)+1) : edim(i);
            pos = [0.1 + params.covariate_gap*(i-1)+sum(w(1:i-1)) 0.1 w(i) 0.85];             
        end   
        
        % main plotting routine
        plot_this_covariate(params,tr{i},bs(:,these_cols),pos,params.covariates_to_plot{i},params.covariate_names{i});
        
        % add cell numbers, group labels and time label if first covariate
        if i==1
            groups = sort(params.group_by);        
            for k=1:params.ngroups
                height = mean(find(groups == k));
                text(-10,height,params.group_labels{k},'color',params.group_colors(k,:),'FontSize',P.font_size,'HorizontalAlignment','right');
            end    
            g=xlabel('Time (s)');  
            if params.show_cell_no
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
 
end

function [bs,edim,tr] = get_data_matrix(stats,params)
    % get combined weight matrix
    for i=numel(params.covariates):-1:1
        [ws,tr{i}] = get_combined_weights_downsample(stats,params.covariates{i},1); % downsample by a factor of 1
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
    params.sig = params.sig(~isnan(params.group_by));

    %% remove unresponsive cells (determined for now by amplitudes of fit coefficients)
    if ~isempty(params.sig)
        bs = bs(params.sig,:);
        params.group_by = params.group_by(params.sig);
    elseif params.responsive_cutoff_prctile>0
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
    set(ax,'clim',params.clim,'ylim',[0.5 params.ncells],P.axes_properties{:});box off;
    ax.YAxis.Visible='off';      
    title(covariate_name,'FontSize',15,'FontWeight','normal');
end

function draw_colorbar(pos,params)
    originalSize = get(gca, 'Position'); % remember current axis position so it can be reset in case colorbar squeezes it        
    h=colorbar('Location','manual');
    ticks = round([0.01:0.01:0.09 0.1:0.1:1.2 2 3 4 5 6 7 8 9 10 20 30 40],2);
    idx = ~ismember(ticks,round([0.01 0.1 0.2 0.5 0.8 0.9 1 1.1 1.2 2 5 10],2));%   [-2 -1 0 1]);
    TickLabels = arrayfun(@num2str,ticks,'uni',0);
    TickLabels(idx) = deal({''});    
    set(h,'FontSize',14,'box','off','LineWidth',0.75,'TickLength',0.05,'Position',[pos(1)+params.covariate_gap+pos(3) 0.66 0.03 0.3],'Ticks',log10(ticks),'TickLabels',TickLabels);
    set(h.Label,'String','Gain','Position',[0.5 0.93 0],'Rotation',0);    
    set(gca,'Position',originalSize);
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