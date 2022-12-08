function gpfa = format_latents(Cells,latents,ref_event,time_window_s,varargin)
    p=inputParser;
    p.KeepUnmatched=true;
    p.addRequired('latents',@(x)validateattributes(x,{'numeric'},{'ndims',3,'nonempty'}));
    p.addRequired('ref_event',@(x)validateattributes(x,{'char'},{'nonempty'}));
    p.addRequired('time_window_s',@(x)validateattributes(x,{'numeric'},{'numel',2,'increasing'}));
    p.addParameter('trial_idx',true(Cells.nTrials,1),@(x)validateattributes(x,{'logical'},{'vector','nonempty'}));
    p.addParameter('nan_trials',false(Cells.nTrials,1),@(x)validateattributes(x,{'logical'},{'vector','nonempty'}));    
    p.addParameter('units',true(numel(Cells.spike_time_s),1),@(x)validateattributes(x,{'logical'},{'vector','nonempty'}));
    p.addParameter('description','',@(x)validateattributes(x,{'char'},{}));
    p.parse(latents,ref_event,time_window_s,varargin{:});
    gpfa = p.Results;
    nan_idx = ismember(find(gpfa.trial_idx),find(gpfa.nan_trials));
    gpfa.score = permute(gpfa.latents,[3 1 2]); % assumes latents are input as Gabe first sent them, as trials x gpfs X times
    gpfa.score(:,nan_idx,:) = NaN;
    gpfa = rmfield(gpfa,'latents');
    gpfa.dim = size(gpfa.score,3);
    gpfa.label = [Cells.rat,' - ',char(Cells.sess_date),' - ',Cells.probe_serial,' - ',num2str(gpfa.dim)',' GPFs'];
    gpfa.repo_path = fileparts(mfilename('fullpath'));
    [gpfa.git_branch,gpfa.commit] = return_git_status(gpfa.repo_path);    
    gpfa.average = 0;
    gpfa.time_s = linspace(gpfa.time_window_s(1),gpfa.time_window_s(2),size(gpfa.score,1));
    gpfa.times = bsxfun(@plus,gpfa.time_s,Cells.Trials.stateTimes.(gpfa.ref_event)(gpfa.trial_idx))';
    gpfa.resolution_s = gpfa.time_s(2) - gpfa.time_s(1);
    gpfa.Trials = Cells.Trials;
    gpfa.kSpikeWindowS = Cells.kSpikeWindowS;
    gpfa.trial_idx = find(gpfa.trial_idx);
end