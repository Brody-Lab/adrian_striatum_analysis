function [cells_mat,params] = MakeDataMatrix(Cells,varargin)
    %
    %MAKEDATAMATRIX make matrix of activity concatenated across trials
    %and cells for data in Cells structure
    %
    %   cells_mat = MAKEDATAMATRIX(Cells) calculates the data matrix using
    %   default parameters and returns it as cells_mat.
    %
    %   MAKEDATAMATRIX(Cells,'ref_event',REF_EVENT) uses data
    %   aligned to event REF_EVENT ('cpoke_in' by default)
    %
    %   MAKEDATAMATRIX(...,'units',UNITS) uses only units specified by
    %   the boolean vector UNITS.
    %
    %   MAKEDATAMATRIX(...,'trial_idx',TRIALS) uses only trials specified by
    %   the boolean vector TRIALS.
    %
    %   MAKEDATAMATRIX(...,'time_edges_s',UNITS) uses binning
    %   specified by the edges given in TIME_EDGES_S in seconds. This corresponds to the
    %   'edges' input of histcounts. (0:0.01:1.5 by default).
    %
    %   cells_mat = MAKEDATAMATRIX(...,'sparse',SPARSE) if SPARSE is true,
    %   cells_mat is returned as a sparse array (useful if memory limited since most values are zero).
    %
    %   cells_mat = MAKEDATAMATRIX(...,'average',AVERAGE) if AVERAGE is true,
    %   cells_mat will be the average firing rate over during ref_event instead of each instance.    
    %
    %   Adrian Bondy, 2021

    %% parse and validate inputs
    p=inputParser;
    p.KeepUnmatched=true;
    ncells = numel(Cells.spike_time_s.cpoke_in);
    p.addParameter('ref_event','cpoke_in',@(x)validateattributes(x,{'char','string'},{'nonempty'}));
    p.addParameter('units',true(ncells,1),@(x)validateattributes(x,{'logical'},{'vector','numel',ncells}));
    p.addParameter('trial_idx',~Cells.Trials.violated,@(x)validateattributes(x,{'logical'},{'vector'}));  
    p.addParameter('resolution_s',0.001,@(x)validateattributes(x,{'numeric'},{'nonnegative'}));
    p.addParameter('time_window_s',[0 1.5],@(x)validateattributes(x,{'numeric'},{'nonnegative'}));    
    p.addParameter('sparse',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('average',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('description','',@(x)validateattributes(x,{'char','string'},{'nonempty'}));    
    p.parse(varargin{:});
    params = p.Results;

    validatestring(params.ref_event,fieldnames(Cells.spike_time_s),'MakeDataMatrix','ref_event');

    params.trial_idx = find(params.trial_idx(:)' & ~isnan(Cells.Trials.stateTimes.(params.ref_event))');
    
    params.units=params.units(:)';
    time_edges_s = params.time_window_s(1):params.resolution_s:params.time_window_s(2);  
    params.time_s = (time_edges_s(1:end-1) + time_edges_s(2:end))/2;    
    % Because clicks are of different number of "trials" alter n_trials for
    % that input
    if strcmp(params.ref_event, 'right_clicks') || strcmp(params.ref_event, 'left_clicks')
        params.trial_idx = 1:length(Cells.spike_time_s.(params.ref_event){1});
    end
    
    %% preallocate and set up loop over cells and trials
    n_time_bins = numel(params.time_s);
    ntrials = numel(params.trial_idx);
    %creates two mats for the average of correct size
    if params.average
            cells_mat = zeros(n_time_bins,sum(params.units));
            hold_mat = zeros(n_time_bins,sum(params.units),ntrials);
    else
            cells_mat = zeros(n_time_bins*ntrials,sum(params.units));
    end
    cell_count=0;    
    
    ref_event_times = Cells.Trials.stateTimes.(params.ref_event)(params.trial_idx);
    params.times = bsxfun(@plus,params.time_s',ref_event_times'); % times of the matrix rows in seconds according to the bcontrol clock
    params.times=params.times(:);
    
    %% loop over cells and trials
    for i=find(params.units)
        spikes_i = Cells.spike_time_s.(params.ref_event){i};
        cell_count = cell_count+1;
        trial_count=0;                
        for j = params.trial_idx
            trial_count = trial_count+1;
            trial_idx  = (1:n_time_bins) + n_time_bins*(trial_count-1);
            if params.average
                hold_mat(:,cell_count,trial_count) = matlab.internal.math.histcounts(spikes_i{j},time_edges_s); 
            else
            % using the hidden internal matlab function instead of "histcounts" avoids repeated argument validation
                cells_mat(trial_idx,cell_count) = matlab.internal.math.histcounts(spikes_i{j},time_edges_s);                     
            end
        end       
    end
    if params.average
        cells_mat = mean(hold_mat,3);
    end   
    %% convert to sparse array if desired
    if params.sparse
        cells_mat = sparse(cells_mat);
    end
    
end
   