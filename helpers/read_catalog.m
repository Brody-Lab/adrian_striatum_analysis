function T = read_catalog(P)
    if nargin<1
        P=get_parameters();
    end
    fprintf('Loading glmfit catalog at %s ...',P.glmfit_catalog_path); tic;   
    opts = detectImportOptions(P.glmfit_catalog_path);
    opts=setvaropts(opts,'save_time','InputFormat','dd-MM-yyyy HH:mm:ss');
    T = readtable(P.glmfit_catalog_path,opts);  
    fprintf(' took %s.\n',timestr(toc));
end