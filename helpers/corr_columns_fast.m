function coef  = corr_columns_fast(x)
    coef = x' * x; % 1/(n-1) doesn't matter, renormalizing anyway
    d = sqrt(diag(coef)); % sqrt first to avoid under/overflow
    coef = coef./d; coef = coef./d'; % coef = coef ./ d*d';
    t = find(abs(coef) > 1); coef(t) = coef(t)./abs(coef(t)); % preserves NaNs
    coef(1:size(x,2)+1:end) = sign(diag(coef)); % preserves NaNs
end
