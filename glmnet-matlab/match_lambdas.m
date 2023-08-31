function fits = match_lambdas(fits)

lambda_intersect = fits{1}.lambda;
for i=1:(numel(fits)-1)
    lambda_intersect = intersect(lambda_intersect,fits{i+1}.lambda);
end

for i=1:numel(fits)
    id = ismember(fits{i}.lambda,lambda_intersect);
    fits{i}.a0 = fits{i}.a0(id);
    fits{i}.df = fits{i}.df(id);
    fits{i}.lambda = fits{i}.lambda(id);
    fits{i}.beta = fits{i}.beta(:,id);
    fits{i}.dim(2) = numel(lambda_intersect);
    fits{i}.dev = fits{i}.dev(id);
end

end