function gpfa = format_latents(Cells,latents,ref_event,time_window_s,varargin)
    p=inputParser;
    p.addRequired('latents',@(x)validateattributes(x,{'numeric'},{'ndims',3,'nonempty'}));
    p.addRequired('ref_event',@(x)validateattributes(x,{'char'},{'nonempty'}));
    p.addRequired('time_window_s',@(x)validateattributes(x,{'numeric'},{'numel',2,'increasing'}));
    p.addParameter('trial_idx',true(Cells.nTrials,1),@(x)validateattributes(x,{'logical'},{'vector','nonempty'}));
    p.addParameter('units',true(numel(Cells.spike_time_s),1),@(x)validateattributes(x,{'logical'},{'vector','nonempty'}));
    p.addParameter('description','',@(x)validateattributes(x,{'char'},{}));
    p.addParameter('make_psth',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.parse(latents,ref_event,time_window_s,varargin{:});
    gpfa = p.Results;
        
    gpfa.score = permute(gpfa.latents,[3 1 2]); % assumes latents are input as Gabe first sent them, as trials x gpfs X times
    gpfa = rmfield(gpfa,'latents');
    gpfa.dim = size(gpfa.score,3);
    gpfa.label = [Cells.rat,' - ',Cells.sess_date,' - ',Cells.probe_serial,' - ',num2str(gpfa.dim)',' GPFs'];
    gpfa.repo_path = fileparts(mfilename('fullpath'));
    [gpfa.git_branch,gpfa.commit] = return_git_status(gpfa.repo_path);    
    gpfa.average = 0;
    gpfa.time_s = linspace(gpfa.time_window_s(1),gpfa.time_window_s(2),1+size(gpfa.score,1));
    gpfa.time_s = (gpfa.time_s(1:end-1) + gpfa.time_s(2:end))/2;
    gpfa.times = bsxfun(@plus,gpfa.time_s,Cells.Trials.stateTimes.(gpfa.ref_event)(gpfa.trial_idx))';
    gpfa.resolution_s = gpfa.time_s(2) - gpfa.time_s(1);
    
    if gpfa.make_psth
        gpfa.psth = get_pc_psth(Cells,gpfa,varargin{:});
    end    
end