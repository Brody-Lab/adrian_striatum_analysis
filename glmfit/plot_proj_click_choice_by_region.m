%% script for breaking down fits by ap group - 3/1/22


P = get_parameters();
cells_table = select_cells('is_in_dorsal_striatum',true);

%%
[glmfit_log,run] = select_glmfit_runs('choice_time_back_s',1.25,'use_trial_history',false,'lambda',30,'include_mono_clicks',true);
[fits,is_responsive] = get_glm_fits(cells_table.recording_name,cells_table.cellno,glmfit_log);
fits_table = [cells_table(is_responsive,:) fits];

[~,idx] = sort(fits_table.AP);
fits_table = fits_table(idx,:);

nboots=500;

for k=1:nboots
    group{k} = equalize_groups(double(fits_table.ap_group),300);        
end
ws_choice = [fits_table.stats.ws];

%%
[glmfit_log,run] = select_glmfit_runs('choice_time_back_s',0,'use_trial_history',false,'lambda',30,'include_mono_clicks',true);
[fits,is_responsive] = get_glm_fits(cells_table.recording_name,cells_table.cellno,glmfit_log);
fits_table = [cells_table(is_responsive,:) fits];

[~,idx] = sort(fits_table.AP);
fits_table = fits_table(idx,:);


ws = [fits_table.stats.ws];

choice_axis = get_choice_axis_from_fits(ws_choice);  
click_axis = get_click_diff_axis_from_fits(ws);


%%
clear h
figure('color',[1 1 1]);
for i=1:numel(P.ap_groups) 
    clear scalar_proj cos_sim choice_norm click_norm
    for k=1:nboots
        [scalar_proj(k,:),scalar_proj_orth(k,:),cos_sim(k,:),choice_norm(k,:),click_norm(k,:)] = get_projection_stats(choice_axis(group{k}==i,:),click_axis(group{k}==i,:));
        
    end
    i
    proj_data(i).scalar_proj = scalar_proj;
    proj_data(i).scalar_proj_orth = scalar_proj_orth;
    proj_data(i).cos_sim = cos_sim;
    proj_data(i).choice_norm = choice_norm;
    proj_data(i).click_norm = click_norm;
    hold on;h(i)=shadedErrorBar(buildGLM.get_tr(ws(1).left_clicks.tr),proj_data(i).scalar_proj,{@mean,@std},{'color',P.ap_group_colors(i,:)});
end
legend([h.mainLine],P.ap_group_labels);
yl=get(gca,'ylim');
set(gca,'ylim',[0 yl(2)],P.axes_properties{:});
ylabel({'Projection of Click Difference (L-R)','Onto Choice Axis'});
xlabel('Time (s)');

clear h
figure('color',[1 1 1]);
for i=1:numel(P.ap_groups)
    hold on;h(i)=plot(mean(proj_data(i).scalar_proj(:,1:25)),mean(proj_data(i).scalar_proj_orth(:,1:25)));
    h(i).Color = P.ap_group_colors(i,:);
    h(i).LineWidth=2;
end
legend(h,P.ap_group_labels);
yl=get(gca,'ylim');
set(gca,'ylim',[0 yl(2)],P.axes_properties{:});
xlabel('Choice Axis Projection');
ylabel('Stimulus Axis Projection');

%% calculate the choice axis trajectory over time, in 3 pcs
clear c choice_axis
for i=1:numel(P.ap_groups)
    parfor k=1:nboots    
        tmp = get_encoding_axis_from_fits(ws(group{k}==i),{'cpoke_out_left','cpoke_out_right'},'time_average',false,'unit_norm',false); 
        [~,c(:,:,k)] = pca(tmp','NumComponents',4);
    end
    choice_axis{i} = c;
    i
end    