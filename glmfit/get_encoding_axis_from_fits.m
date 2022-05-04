function choice_axis = get_encoding_axis_from_fits(stats,covariate,varargin)
    p=inputParser;
    p.addRequired('covariate',@(x)validateattributes(x,{'char','cell'},{'nonempty'}));
    p.addParameter('time_lim_s',[0 1]);
    p.addParameter('time_average',false);
    p.addParameter('unit_norm',false);    
    p.parse(covariate,varargin{:});
    params = p.Results;
    
    if iscell(covariate)
        if numel(covariate)==1
            [choice_axis,tr] = get_combined_weights_downsample(stats,covariate{1},10);            
        elseif numel(covariate)==2          % if coviate is a two-elements cell array, take the difference (useful for left-right difference encoding)   
            [ws_left,tr] = get_combined_weights_downsample(stats,covariate{1},10);
            ws_right = get_combined_weights_downsample(stats,covariate{2},10);
            choice_axis = ws_right.data - ws_left.data;
        else
            error('Cell array input for "covariate" must have size 1 or 2.');
        end
    else
        [choice_axis,tr] = get_combined_weights_downsample(stats,covariate,10);                    
    end
    
    idx = tr > params.time_lim_s(1) & tr < params.time_lim_s(2);
    
    if params.time_average
        choice_axis = mean(choice_axis(:,idx),2);
    end

    if params.unit_norm
        n = sqrt(sum(choice_axis.^2,1));
        choice_axis = bsxfun(@times,choice_axis,1./n);
    end
   
end