%% script for breaking down fits by ap group - 3/1/22


P = get_parameters();
[glmfit_log,run] = select_glmfit_runs('choice_time_back_s',0.5,'use_trial_history',false,'lambda',30,'include_mono_clicks',true);
cells_table = select_cells('is_in_dorsal_striatum',true);
[fits,is_responsive] = get_glm_fits(cells_table.recording_name,cells_table.cellno,glmfit_log);
fits_table = [cells_table(is_responsive,:) fits];

[~,idx] = sort(fits_table.AP);
fits_table = fits_table(idx,:);

nboots=500;

for k=1:nboots
    group{k} = equalize_groups(double(fits_table.ap_group),300);        
end
ws = [fits_table.stats.ws];

figure;
for i=1:numel(P.ap_groups) 
    clear scalar_proj cos_sim choice_norm click_norm
    parfor k=1:nboots
        [scalar_proj(k,:),cos_sim(k,:),choice_norm(k,:),click_norm(k,:)] = get_proj_click_choice(ws(group{k}==i));
        
    end
    i
    proj_data(i).scalar_proj = scalar_proj;
    proj_data(i).cos_sim = cos_sim;
    proj_data(i).choice_norm = choice_norm;
    proj_data(i).click_norm = click_norm;
    hold on;shadedErrorBar(buildGLM.get_tr(ws(1).left_clicks.tr),proj_data(i).scalar_proj,{@mean,@std},{'color',P.ap_group_colors(i,:)});
end

