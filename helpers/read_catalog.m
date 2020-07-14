function T = read_catalog()
    P=get_parameters();
    opts = detectImportOptions(P.glmfit_catalog_path);
    opts=setvaropts(opts,'save_time','InputFormat','dd-MM-yyyy HH:mm:ss');
    T = readtable(P.glmfit_catalog_path,opts);  
end