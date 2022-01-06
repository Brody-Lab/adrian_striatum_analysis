function [keep_idx,n_template,n_matched,bin_edges] = equalize_trials(features,template_idx,nbins)

    % features is a matrix of size ntrials x nfeatures describing the
    % values of each trial in the feature space to be equalized over. There
    % is no limit to the dimensionality of this array.
    
    % template_idx is the boolean vector with 1s corresponding to trials
    % that are in the group whose distributed is to be used as the
    % template.
    
    % nbins is a vector of length nfeatures describing the number of
    % bins per feature
    
    % output keep_idx indexes the trials in the matching set to keep to
    % equalize conditions with the template set
    
    
    %% parse and validate inputs
    p=inputParser;
    p.addRequired('features',@(x)validateattributes(x,{'numeric'},{'ndims',2,'nonempty'}));
    p.addRequired('template_idx',@(x)validateattributes(x,{'logical'},{'vector','nonempty'}));    
    p.addRequired('nbins',@(x)validateattributes(x,{'numeric'},{'vector','nonempty'}));
    p.parse(features,template_idx,nbins);
    
    [ntrials,nfeatures] = size(features);
    if numel(template_idx)~=ntrials
        error('Number of elements of template_idx must match size of first dimension of features');
    end
    if numel(nbins)~=nfeatures
        error('Number of elements of nbins must match size of second dimension of features');
    end
    
    %% get bin edges and discretize into bins
    for i=nfeatures:-1:1
        bin_edges{i} = linspace(min(features(template_idx,i)),max(features(template_idx,i)),nbins(i)+1);
        features_discretized(:,i) = discretize(features(:,i),bin_edges{i}); 
    end
    
    %% find all present combinations of features
    nans = any(isnan(features_discretized),2);
    u = unique(features_discretized(~nans,:),'rows');
    
    %% for each combination, find number of template matches and sample the same number in the matching set
    if nfeatures>1  
        [n_template,n_matched] = deal(zeros(nbins));
    else
        [n_template,n_matched] = deal(zeros(nbins,1));        
    end
    keep_idx=[];
    for i=1:size(u,1)
        if nfeatures>1
            tmp=num2cell(u(i,:));            
            idx = sub2ind(nbins,tmp{:});
        else
            idx=u(i);
        end 
        n_template(idx) = sum(ismember(features_discretized(template_idx,:),u(i,:),'rows'));
        matched = find(~template_idx & ismember(features_discretized,u(i,:),'rows'));
        n_matched(idx) = numel(matched);
        if n_matched(idx)<n_template(idx)
            warning('Only %g of %g matched for bin %g.\n',n_matched(idx),n_template(idx),i);
            nsamples = n_matched(idx);
        else
            nsamples = n_template(idx);
        end
        keep_idx = union(keep_idx,matched( randsample(n_matched(idx),nsamples)));
    end           
end