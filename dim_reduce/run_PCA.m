function stats = run_PCA(tot_cell_mat,varargin) %dim,norm_factor)
    % normalizes created unit by time sample and runs PCA, returning PCs,
    % scores, and percentage explained, and the normalized activity mat.
    %
    %[coeff,score,explained,norm_cell_mat] = run_PCA(tot_cell_mat)
    % calculates PCs for neural data and returns key metrics and normalized
    % data
    %   
    %run_PCA(tot_cell_mat,'dim',DIM) indicates which dim space should be
    %used to calculate PCs 1 is neural 2 is time
    %
    %run_PCA(tot_cell_mat,'norm_factor',NORM) specifies a value for soft
    %normalization in the denominator 
    
    p=inputParser;
    p.KeepUnmatched=true;
    p.addRequired('tot_cell_mat',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
    p.addParameter('norm_factor',5,@(x)validateattributes(x,{'numeric'},{'nonnegative'}));
    p.addParameter('npcs',Inf,@(x)validateattributes(x,{'numeric'},{'nonnegative'}));    
    p.parse(tot_cell_mat,varargin{:});
    params = p.Results; 

    %current normalization subtract columns then divide by range
    norm_cell_mat = tot_cell_mat -mean(tot_cell_mat);
    norm_cell_mat = norm_cell_mat./(range(norm_cell_mat)+params.norm_factor);
    
    if ismember('npcs',p.UsingDefaults) 
        params.npcs = size(norm_cell_mat,2);
    elseif size(norm_cell_mat,2)<params.npcs
        error('You requested more PCs than colums of X.');
    end
    
    %save coeff, weightings, and var explained
    [stats.coeff, stats.score,stats.latent, stats.tsquared, stats.explained] = pca(norm_cell_mat,'NumComponents',params.npcs);
    
end
