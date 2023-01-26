% Matlab script for regenerating the processed data for problem set 1 for NEU 437 Spring 2023.
% this is a dataset of rat ADS neurons recorded by Adrian Bondy
% during performance of the Poisson Clicks task using a Neuropixels 1.0 probe.
% for more information about this data, visit: https://github.com/Brody-Lab/adrian_striatum_analysis

% Note: To run this script, it is necessary to have the Brody Lab repos
% "adrian_striatum_analysis" and "labwide_pbups_analysis" in your matlab
% path.

% AGB 1/26/2023

save_path = "pset1_data";
recording_name = "A249_2020_09_07";
Cells = load_Cells_file(recording_name);
exclude_trials = validate_trials(Cells.Trials,'mode','agb_glm');
params = get_pcs(Cells,'resolution_s',5e-3,'trial_idx',~exclude_trials,...
    'exclude_cells',~Cells.is_in_dorsal_striatum);
fields=fieldnames(params);
mask = {'ref_event','resolution_s','time_window_s',...
    'trial_idx','units','times','cells_mat','sessid','rat','sess_date'}; % fields to keep
params = rmfield(params,fields(~ismember(fields,mask)));
save(save_path,'-struct','params');


