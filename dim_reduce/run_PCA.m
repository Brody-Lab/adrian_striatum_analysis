function [pca_output,params] = run_PCA(tot_cell_mat,varargin) %dim,norm_factor)
    % normalizes created unit by time sample and runs PCA, returning PCs,
    % scores, and percentage explained, and the normalized activity mat.
    %
    %[coeff,score,explained,norm_cell_mat] = run_PCA(tot_cell_mat)
    % calculates PCs for neural data and returns key metrics and normalized
    % data
    %
    %run_PCA(tot_cell_mat,'norm_factor',NORM) specifies a value for soft
    %normalization in the denominator 
    %
    %run_PCA(tot_cell_mat,'npcs',NPCS) specifies the number of principal
    %components to use
    %
    %run_PCA(tot_cell_mat,'nsubsamples',NSUBSAMPLES) specifies the number of
    %times to subsample the population of cells (the size of the subsampled
    %population is determined by NPCS
    %
    % N.B. Can go faster if tot_cell_mat is a gpuArray and/or single
    % precision.
    
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('norm_factor',5,@(x)validateattributes(x,{'numeric'},{'nonnegative'}));
    p.addParameter('npcs',Inf,@(x)validateattributes(x,{'numeric'},{'positive'}));  
    p.addParameter('nsubsamples',1,@(x)validateattributes(x,{'numeric'},{'positive'}));        
    p.parse(varargin{:});
    params = p.Results; 

    %current normalization subtract columns then divide by range
    norm_cell_mat = tot_cell_mat -mean(tot_cell_mat);
    norm_cell_mat = norm_cell_mat./(range(norm_cell_mat)+params.norm_factor);
    params.ncells = size(norm_cell_mat,2);
    
    if ismember('npcs',p.UsingDefaults) 
        params.npcs = params.ncells;
    elseif params.ncells<params.npcs
        error('You requested more PCs than colums of X.');
    end
    
    if params.nsubsamples>1
        if params.npcs==params.ncells
            error('To set nsubsamples>1, npcs must be smaller than the number of cells (i.e. size of tot_cell_mat in dimension 2).');
        end
        for i=params.nsubsamples:-1:1
            pca_output(i).sidx = randsample(params.ncells,params.npcs);
            [pca_output(i).coeff, pca_output(i).score,pca_output(i).latent] = ...
                pca(norm_cell_mat(:,pca_output(i).sidx),'NumComponents',params.npcs);        
            pca_output(i).explained = 100*pca_output(i).latent/sum(pca_output(i).latent);
        end
    else
        [pca_output.coeff, pca_output.score,pca_output.latent] = pca(norm_cell_mat,'NumComponents',params.npcs);
        pca_output.explained = 100*pca_output.latent/sum(pca_output.latent);        
    end
end
