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
    p.addParameter('kPETH',get_PETH_params('resolution_s',0.05,'type','GAUSS'));
    % you can pass mask_states directly to get_psth_from_Cells
    p.parse(varargin{:});
    params=p.Results;
    if isempty(params.column_name)
        params.column_name_root = [params.ref_event,'_peth'];
    else
        params.column_name_root = params.column_name;   
    end
    params.column_name = {[params.column_name_root,'_observed'] [params.column_name_root,'_predicted'] [params.column_name_root,'_observed_se'] [params.column_name_root,'_predicted_se']};
    fprintf('Using "%s" and "%s" as the column names.\n',params.column_name{1},params.column_name{2});      
    recording_names = unique(fits_table.recording_name);
    run = unique(fits_table.run);
    if ~isscalar(run)
        error('Currently we assume a single run for the whole fits table for convenience although this is not a necessary assumption here.');
    end
    P=get_parameters();
    fit_path=P.fit_path;
    
    for i=numel(recording_names):-1:1
        these_rows{i} = fits_table.recording_name == recording_names(i);
        cellnos{i} = fits_table.cellno(these_rows{i});
        stats{i} = fits_table.stats(these_rows{i});   
    end
    
    parfor i=1:numel(recording_names)
        params_file = load(fullfile(fit_path,'fits',recording_names(i),run,'glmfit_params.mat'));        
        Cells = load_Cells_file(recording_names(i));
        [~,p,spikes] = fit_glm_to_Cells(Cells,params_file.params,'fit',false,'keepX',true,'save',false,'cellno',cellnos{i});        
        [~,idx] = ismember(cellnos{i},[spikes.cellno]);
        spikes=spikes(idx);
        % predict Y_hat using X in params and stats for each cell
        for c=1:numel(spikes)
            if isfield(stats{i}(c),'B')
                spikes(c).Yhat = p.link.Inverse([p.X p.dm{c}.X]*stats{i}(c).B(2:end) + stats{i}(c).B(1));
            else
                spikes(c).Yhat = p.link.Inverse([p.X p.dm{c}.X]*stats{i}(c).beta(2:end) + stats{i}(c).beta(1));                
            end
        end    
        [peth_observed{i},peth_predicted{i},peth_observed_se{i},peth_predicted_se{i},boot_sim_observed{i},boot_sim_predicted{i}] = get_psth_from_glmfit(p.dm_base.dspec.expt,spikes,params.ref_event,'kPETH',params.kPETH,'separate_by',params.separate_by);    
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
    end

end





