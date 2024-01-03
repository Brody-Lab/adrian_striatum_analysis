function [fits_table,kPETH] = add_peth_to_fits_table(fits_table,varargin)


    % similar to add_peth_to_cells_table EXCEPT it doesn't use spike times from Cells files,
    % it uses binned spikes as used in fit_glm_to_Cells. Makes for clean
    % comparison of predicted and observed PETHs under glm model.

    % 
    
    %% parse inputs
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('ref_event','cpoke_in',@(x)validateattributes(x,{'char'},{'nonempty'}));    
    p.addParameter('separate_by','signal_strength',@(x)validateattributes(x,{'char'},{}));
    p.addParameter('column_name','',@(x)validateattributes(x,{'char'},{}));    
    p.addParameter('kPETH',get_PETH_params('std_s',0.1,'std_s_clicks',0.1));
    p.addParameter('nresamples',100); % number of bootstrap resamples    
    % you can pass mask_states directly to get_psth_from_Cells
    p.parse(varargin{:});
    params=p.Results;
    kPETH=params.kPETH;
    if isempty(params.column_name)
        params.column_name_root = [params.ref_event,'_peth'];
    else
        params.column_name_root = params.column_name;   
    end
    params.column_name = {[params.column_name_root,'_observed'] [params.column_name_root,'_predicted'] [params.column_name_root,'_observed_se'] [params.column_name_root,'_predicted_se']};
    fprintf('Using "%s" and "%s" as the column names.\n',params.column_name{1},params.column_name{2});      
    recording_names =   unique(fits_table.recording_name);%(fits_table.sess_date==datetime("20-Jul-2023")));

    run = unique(fits_table.run);
    if ~isscalar(run)
        warning('More than one run in the fits table. If spikes are generated in different ways across these runs, peth generation will be affected, i.e. time window of trial or bin_size_s');
        run=run(1);
    end
    P=get_parameters();
    fit_path=P.fit_path;
    
    for i=numel(recording_names):-1:1
        these_rows{i} = fits_table.recording_name == recording_names(i);
        cellnos{i} = fits_table.cellno(these_rows{i});
        stats{i} = fits_table.stats(these_rows{i});   
    end
    for i=1:numel(recording_names)
        params_file = load(fullfile(fit_path,'fits',recording_names(i),run,'glmfit_params.mat'));        
        Cells = load_Cells_file(recording_names(i));
        Cells=add_first_click_state(Cells);        
        exclude_trials = validate_trials(Cells.Trials,'mode','agb_glm','quiet',true,'require_clicks',true);
        [choice_MI{i},choice_MI_pval{i}] = get_pref_choice(Cells,cellnos{i},'exclude_trials',exclude_trials);        
        [~,p,spikes] = fit_glm_to_Cells(Cells,params_file.params,'fit',false,'keepX',true,'save',false,'cellno',cellnos{i});        
        [~,idx] = ismember(cellnos{i},[spikes.cellno]);
        spikes=spikes(idx);
        % predict Y_hat using X in params and stats for each cell
        s=stats{i};
        for c=1:numel(spikes)
            if isfield(s(c),'B')
                spikes(c).Yhat = p.link.Inverse([p.X p.dm{c}.X]*s(c).B(2:end) + s(c).B(1));
            else
                spikes(c).Yhat = p.link.Inverse([p.X p.dm{c}.X]*s(c).beta(2:end) + s(c).beta(1));                
            end
        end    
        peth_mean_rate{i} = p.totalSpikes ./ size(p.X,1) ./ p.bin_size_s;
        [peth_observed{i},peth_predicted{i},peth_observed_se{i},peth_predicted_se{i},boot_sim_observed{i},~] = get_psth_from_glmfit(p.dm_base.dspec.expt,spikes,params.ref_event,varargin{:},'kPETH',params.kPETH,'separate_by',params.separate_by);    
        for c=1:numel(cellnos{i})
            r2{i}(c) = rsquare(peth_observed{i}{c}(:),peth_predicted{i}{c}(:));
            r{i}(c) = corr(peth_observed{i}{c}(:),peth_predicted{i}{c}(:));
        end            
    end
    
    for i=1:numel(recording_names)
        
        fits_table.(params.column_name{1})(these_rows{i}) = peth_observed{i};
        fits_table.(params.column_name{2})(these_rows{i}) = peth_predicted{i};       
        fits_table.(params.column_name{3})(these_rows{i}) = peth_observed_se{i};
        fits_table.(params.column_name{4})(these_rows{i}) = peth_predicted_se{i};            
        
        fits_table.([params.column_name_root,'_r2'])(these_rows{i}) = r2{i};
        fits_table.([params.column_name_root,'_r'])(these_rows{i}) = r{i};    
        fits_table.([params.column_name_root,'_boot_sim_observed'])(these_rows{i}) = boot_sim_observed{i};     
        fits_table.([params.column_name_root,'_mean_rate'])(these_rows{i}) = peth_mean_rate{i};
        fits_table.choice_MI(these_rows{i}) = choice_MI{i};
        fits_table.choice_MI_pval(these_rows{i}) = choice_MI_pval{i};             
    end

end





