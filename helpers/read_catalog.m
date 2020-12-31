function T = read_catalog(P)
    if nargin<1
        P=get_parameters();
    end
    fprintf('Loading glmfit catalog at %s ...',P.glmfit_catalog_path); tic;   
    T = readtable(P.glmfit_catalog_path);  
    fprintf(' took %s.\n',timestr(toc));
end