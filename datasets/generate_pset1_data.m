% Matlab script for regenerating the processed data for problem set 1 for NEU 437 Spring 2023.
% this is a dataset of rat ADS neurons recorded by Adrian Bondy
% during performance of the Poisson Clicks task using a Neuropixels 1.0 probe.
% for more information about this data, visit: https://github.com/Brody-Lab/adrian_striatum_analysis

% Note: To run this script, it is necessary to have the Brody Lab repos
% "adrian_striatum_analysis" and "labwide_pbups_analysis" in your matlab
% path.

% AGB 1/26/2023

save_path = "pset1_data";
recording_name = "T219_2019_12_20";
Cells = load_Cells_file(recording_name);
Cells = add_first_click_state(Cells);
exclude_trials = validate_trials(Cells.Trials,'mode','agb_glm');
clicks_on = Cells.Trials.stateTimes.first_click - Cells.Trials.stateTimes.cpoke_in;
clicks_on = clicks_on(~exclude_trials);
params = get_pcs(Cells,'resolution_s',2e-2,'trial_idx',~exclude_trials,...
    'exclude_cells',~Cells.is_in_dorsal_striatum);
params.clicks_on=round(0.5+clicks_on/params.resolution_s);
params.gamma = Cells.Trials.gamma(~exclude_trials);
fields=fieldnames(params);
params.rat=char(params.rat);
mask = {'ref_event','resolution_s','time_window_s',...
    'trial_idx','units','time_s','cells_mat','sessid','rat','sess_date','times','clicks_on','gamma'}; % fields to keep
params = rmfield(params,fields(~ismember(fields,mask)));
params.spikes = reshape(params.cells_mat,[size(params.times) size(params.cells_mat,2)]);
params = rmfield(params,{'cells_mat','times'});
params.sess_date=char(params.sess_date);
params.spikes=permute(params.spikes,[2 1 3]); % now spikes is ntrials X ntimebins X ncells
params.n_left_clicks = cellfun(@numel,Cells.Trials.leftBups(~exclude_trials));
params.n_right_clicks = cellfun(@numel,Cells.Trials.rightBups(~exclude_trials));
params.went_right = Cells.Trials.pokedR(~exclude_trials);
save(save_path,'-struct','params');

