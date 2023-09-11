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
[~,P.hostname] = system('hostname');
P.hostname=deblank(P.hostname);
P.pc_data_path = fullfile('D:','adrian_striatum_analysis');
P.jukebox_data_path = fullfile('X:','abondy','adrian_striatum_analysis');
P.tiger_data_path = '/scratch/gpfs/abondy/adrian_striatum_analysis';
P.della_data_path = '/scratch/gpfs/abondy/adrian_striatum_analysis';
P.spock_data_path = '/mnt/cup/labs/brody/abondy/adrian_striatum_analysis';
P.on_cluster=true;
if strncmp(P.hostname,'tiger',5)
    P.data_path = P.tiger_data_path;
    P.fit_path = P.tiger_data_path;      
elseif strncmp(P.hostname,'della',5)
    P.data_path = P.della_data_path;
    P.fit_path = P.della_data_path;  
elseif strncmp(P.hostname,'spock',5) || strncmp(P.hostname,'scotty',5)
    P.data_path = P.spock_data_path;
    P.fit_path = P.spock_data_path;       
else
    P.on_cluster=false;
    P.data_path = P.pc_data_path;   
    P.fit_path = P.pc_data_path;
end
if ~isfolder(P.data_path)
   error('Data path does not exist: %s',P.data_path); 
end
P.recordings_path = fullfile(P.repository_path,'recordings_log.csv');
P.cells_table_path = fullfile(P.data_path,'cells_table.mat');
P.sessions_table_path = fullfile(P.repository_path,'sessions_table.mat');
P.histology_table_path = fullfile(P.repository_path,'histology_table.csv');
P.implant_log_path = fullfile(P.repository_path,'implant_log.json');


%% analysis
P.ap_groups = {[1.5 Inf],[0 1.5],[-1.5 0],[-Inf -1.5]};

%% glm fitting
P.glmfit_catalog_path = fullfile(P.data_path,'glmfit_log.mat');
% if all these parameters are the same for a cell's fit, the fits should be
% identical, unless a code change produced changes in the fitting algorithm
P.glmfit_catalog_keys = {'rat','recording_name','sess_date','sessid','bin_size_s','fit_method','relaxed','alpha','lambda','kfold','distribution','link','choice_time_back_s',...
    'git_branch','git_commit','include_stereo_click','phi','tau_phi','fit_adaptation','within_stream',...
    'include_mono_clicks','dm_scaling_mode','use_trial_history','isi_cutoff','maxIter','nClickBins','separate_clicks_by',...
    'confine_clicks_to_stim_period','confine_coupling_to_stim_period','coupling_cells','stim_period_boundary','stim_period_mode','separate_clicks_by_outcome','lambda_correct_zscore',...
    'responsiveness_boundary','responsiveness_boundary_offset','ref_event','separate_cpoke_out_by_choice','separate_cpoke_in_by_choice','boundary_mode','clicks_boundary','click_boundary_offset',...
    'coupling_boundary','coupling_boundary_offset','foldseed'}; 
P.glmfit_catalog_params = [P.glmfit_catalog_keys 'run','hostname','save_time','runtime_sec','saved_cells','responsive_cells',...
    'n_missing_cells','minSpkParamRatio','minResponsiveFrac','dm','glmnet_thresh','max_glmnet_attempts','glmnetOptions'];

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
P.ap_group_labels_xtick = ["     site 1\newline      ADS" "site 2" "site 3" "     site 4\newline        TS"];
P.ap_group_labels = ["site 1 (ADS)" "site 2" "site 3" "site 4 (TS)"];
P.pref_color = [0 0 0];
P.nonpref_color = [1 0 0];
P.covariate_order = {'previous_left_hit','previous_right_hit','previous_left_error','previous_right_error',...
                     'cpoke_in','stereo_click','left_clicks','right_clicks','stereo_click_hit','stereo_click_error','left_clicks_hit','left_clicks_error',...
                     'right_clicks_hit','right_clicks_error',...
                     'left_clicks1','left_clicks2','left_clicks3',...
                     'right_clicks1','right_clicks2','right_clicks3',...
                     'cpoke_out_left','cpoke_out_right','spoke_left_hit','spoke_right_hit',...
                     'spoke_left_miss','spoke_right_miss'};
P.covariate_names = {'Previous Left Correct','Previous Right Correct','Previous Left Error','Previous Right Error',...
    'Initiate Fixation','First (Stereo) Click','Left Clicks','Right Clicks','First (Stereo) Click Correct','First (Stereo) Click Error','Left Clicks Correct','Left Clicks Error','Right Clicks Correct','Right Clicks Error',...
    'Left Clicks 1','Left Clicks 2','Left Clicks 3','Right Clicks 1','Right Clicks 2','Right Clicks 3','Left Choice','Right Choice','Left Correct',...
    'Right Correct','Left Error','Right Error'};                 
%P.covariate_group_id = [1 1 1 1 2 3 3 3 4 4 5 5 6 6]; not used or updated
%with current set of names
%P.covariate_group_names = {'Trial History','Initiate Fixation','Clicks','Choice','Hit','Miss'};
P.gamma_color_groups = flip(copper(4));
P.gamma_ranges = [-5 -2 0 2 5];
P.lambdas = logspace(-1,2,4);
P.signal_strength_percentiles = 0:25:100;
P.gamma_labels={'strong preferred','weak preferred','weak null','strong null'};

P.example_cell_recording_names = ["T219_2019_12_09"  "A243_2019_10_29" "T219_2019_12_05" "A249_2020_09_11" "A256_2020_09_23" "A249_2020_09_25"  ];
P.example_cellno = [12 111 71 188 351 168];

               