% GET_PARAMETERS Returns a structure with the parameters for analyses, such as file paths
%
%=OUTPUT
%
%   P
%       A structure with parameters for analyses
%
function P = get_parameters()
%% paths
P.repository_path = fileparts(mfilename('fullpath'));
P.recordings_path = fullfile(P.repository_path,'recordings_log.csv');
[~,P.hostname] = system('hostname');
P.hostname=deblank(P.hostname);
P.pc_data_path = fullfile('D:','adrian_striatum_analysis');
P.jukebox_data_path = fullfile('X:','abondy','adrian_striatum_analysis');
P.tiger_data_path = '/scratch/gpfs/abondy/adrian_striatum_analysis';
if strncmp(P.hostname,'tiger',5)
    P.on_tiger=true;
    P.data_path = P.tiger_data_path;
    P.fit_path = P.tiger_data_path;    
else
    P.on_tiger=false;
    P.data_path = P.pc_data_path;   
    P.fit_path = P.jukebox_data_path;
end
if ~isfolder(P.data_path)
   error('Data path does not exist: %s',P.data_path); 
end
P.cells_table_path = fullfile(P.data_path,'cells_table.mat');
P.sessions_table_path = fullfile(P.data_path,'sessions_table.mat');

%% analysis
P.ap_groups = {[1.5 Inf],[0 1.5],[-1.5 0],[-Inf -1.5]};

%% glm fitting
P.glmfit_catalog_path = fullfile(P.data_path,'glmfit_log.mat');
% if all these parameters are the same for a cell's fit, the fits should be
% identical, unless a code change produced changes in the fitting algorithm
P.glmfit_catalog_keys = {'recording_name','phi','tau_phi','fit_adaptation','bin_size_s','include_mono_clicks','dm_scaling_mode',...
    'choice_time_back_s','kfold','distribution','link','within_stream','lambda','git_branch','git_commit','include_stereo_click','maxIter'}; 
P.glmfit_catalog_params = [P.glmfit_catalog_keys 'rat','sess_date','sessid','run','hostname','save_time','runtime_sec','saved_cells','responsive_cells','n_missing_cells','minSpkParamRatio','minResponsiveFrac','dm'];

%% tagging
P.tagging_metrics={'reliability','dp','tp','signrank','auc','mi'};
P.tagging_metric_names = {'reliability','d-prime','p-value of paired t-test','p-value of sign rank test','area under ROC curve','modulation index'};

%% plotting
P.figure_image_format = {'png'};
P.figure_position= [rand*1000, rand*1000, 1000*1.1, 850*1.1];
P.font_size = 14;
P.axes_properties = {'FontSize', P.font_size, ...
                     'Color', 'none', ...
                     'TickDir', 'Out',...
                     'Nextplot', 'add', ...
                     'LineWidth', 1,...
                     'XColor', [0,0,0], ...
                     'YColor', [0,0,0], ...
                     'LabelFontSizeMultiplier', 1,...
                     'box','off'};  
P.panel_label_font_size = P.font_size * 1.5;
P.panel_label_pos = [0.1, 0.9, 0.1, 0.1];
P.panel_labels = char(65:90);
P.ap_group_colors = copper(4);
P.ap_group_colors = cool(4);
P.ap_group_colors = parula(6);
P.ap_group_colors = P.ap_group_colors(2:5,:);
%P.ap_group_colors = P.ap_group_colors([1 12 18 33],:);
% P.ap_group_colors = redblue(8);
% P.ap_group_colors = P.ap_group_colors([1 4 5 8],:);
P.ap_group_labels_xtick = ["     site 1\newlinemost anterior" "site 2" "site 3" "     site 4\newlinemost posterior"];
P.ap_group_labels = ["site 1 (most anterior)" "site 2" "site 3" "site 4 (most posterior)"];
P.pref_color = [0 0 0];
P.nonpref_color = [1 0 0];
P.covariate_order = {'cpoke_in','stereo_click','left_clicks','right_clicks',...
                     'cpoke_out_left','cpoke_out_right','spoke_left_hit','spoke_right_hit',...
                     'spoke_left_miss','spoke_right_miss'};
P.covariate_names = {'Enter Center Port','First (Stereo) Click','Left Clicks','Right Clicks','Left Choice','Right Choice','Hit Left','Hit Right','Miss Left','Miss Right'};                 
P.covariate_group_id = [1 2 2 2 3 3 4 4 5 5];
P.covariate_group_names = {'Enter Center Port','Clicks','Choice','Hit','Miss'};
P.gamma_color_groups = flip(copper(4));
               