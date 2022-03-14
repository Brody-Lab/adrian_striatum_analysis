function psth = get_psth(Cells,cellno,varargin)
    p = inputParser;
    p.KeepUnmatched=true;
    p.addParameter('trial_idx',~Cells.Trials.violated,@(x)validateattributes(x,{'logical','numeric'},{'vector'}));  
    p.addParameter('nresamples',1); 
    p.addParameter('states',{},@(x)validateattributes(x,{'cell'},{}));
    p.addParameter('kPETH',get_PETH_params());
    p.parse(varargin{:});
    params = p.Results;
    kPETH = params.kPETH; 
    % only include desired trials in psth calculation
    if islogical(params.trial_idx)
        params.trial_idx = find(params.trial_idx);
    end
    
    %% make PSTHs    
    if ismember('states',p.UsingDefaults)
        states=kPETH.refEvents;
    else
        states = params.states;
    end    
    for eve=string(states(:)')  
        for c=1:numel(cellno)
            psth(c).(eve)=calc_psth(Cells.spike_time_s.(eve){cellno(c)}(params.trial_idx),kPETH.timeS.(eve),...
                kPETH.type.(eve),kPETH.stdS.(eve));
            if params.nresamples>1
                psth(c).(eve) = gather(bootstrp(params.nresamples,@nanmean,gpuArray(psth(c).(eve))));              
            end    
        end
    end
end