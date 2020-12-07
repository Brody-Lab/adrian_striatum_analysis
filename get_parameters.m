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
P.pc_data_path = fullfile('X:','abondy','adrian_striatum_analysis');
P.tiger_data_path = '/tigress/abondy';
if strncmp(P.hostname,'tiger',5)
    P.on_tiger=true;
    P.data_path = P.tiger_data_path;
else
    P.on_tiger=false;
    P.data_path = P.pc_data_path;    
end
if ~isdir(P.data_path)
   error('Data path does not exist: %s',P.data_path); 
end
P.glmfit_catalog_path = fullfile(P.data_path,'glmfit_log.csv');
P.glmfit_catalog_keys = {'phi','tau_phi','fit_adaptation','bin_size_s','include_mono_clicks','choice_time_back_s','kfold','distribution','link','within_stream','lambda'}; % if all these parameters are the same for a cell's fit, the fits should be identical, except perhaps for very minor changes in the fitting algorithm.
P.glmfit_catalog_params = [P.glmfit_catalog_keys 'rat','sess_date','sessid','run','cells_file','git_branch','git_commit','hostname','save_time'];
    

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
                     'ActivePositionProperty','Position',...
                     'box','off'};  
P.panel_label_font_size = P.font_size * 1.5;
P.panel_label_pos = [0.1, 0.9, 0.1, 0.1];
P.panel_labels = char(65:90);

               