function avg_corr = cluster_consistency(X,max_clust) 
n_rep=20;
n_samples= size(X,1);
metric = 'euclidean';
method = 'ward';

for i=1:n_rep
    idx= randsample(n_samples,n_samples,true);
    
    %% determine tree structure
    tree = linkage(X(idx,:),method,metric);
    
    %% determine distance cutoff to produce desired number of clusters   
    clusterfun = @(cutoff)cluster(tree,'Cutoff',cutoff,'Criterion','distance'); 
    [~,T(:,i)] = find_cutoff(clusterfun,max_clust,5,10);
    
end
    avg_corr = get_clustering_corr(X,T);
end
    
    
        
        
        
function avg_corr = get_clustering_corr(X,T)
        
n_clust = numel(unique(T));
for i = 1:size(T,2)
    for g= 1:n_clust
        means(i,g,:) = mean(X(T(:,i)==g,:));
    end
end
        
count=0;
for i=1:size(T,2)
    for k=1:size(T,2)
        if k<=i
            continue
        else
            count=count+1;
            corrs=corr(squeeze(means(k,:,:))',squeeze(means(i,:,:))');
            avg_corr(count,:) = (max(corrs));
        end
    end
end        
end


function [cutoff,T] = find_cutoff(clusterfun,max_clust,guess,step)
    last_guess=guess;
    min_guess = eps;
    count=0;
    while true
        T = clusterfun(guess);
        n_clust = numel(unique(T));
        count=count+1;
        if n_clust>max_clust            
            if guess<last_guess
                step=step/2;
            end
            last_guess=guess;            
            guess=guess+step;
        elseif n_clust==max_clust
            cutoff=guess;
            fprintf('Found cutoff of %g for %g clusters in %g iterations.\n',cutoff,n_clust,count);
            return
        elseif n_clust<max_clust
            if guess>last_guess
                step=step/2;
            end
            last_guess=guess;            
            guess=guess-step;
        end
        if guess<min_guess
            guess=eps;
        end        
    end
end
