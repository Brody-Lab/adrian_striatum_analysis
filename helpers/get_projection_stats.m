function [scalar_proj,cos_sim,norma,normb] = get_projection_stats(a,b)
    % project vector B onto A
    % if B is a matrix, project each column separately
    if isempty(a) || isempty(b)
        error('A and B can''t be empty.');
    end
    if size(a,1) ~= size(b,1)
        error('A and B must have the same number of rows, i.e. dimension');
    end
    if size(a,2)~=1
        error('A must be a column vector');
    end
    dotprod = a'*b;
    norma = norm(a);
    normb = sqrt(sum(b.^2));
    cos_sim = dotprod ./ (norma * normb);
    scalar_proj = cos_sim .* normb;
end