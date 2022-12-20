function histology_table = load_histology_table()
    %% imports the histology .csv file into matlab
    P = get_parameters;
    opts = detectImportOptions(P.histology_table_path,'TextType','string');
    opts = opts.setvartype('acquisition_date','datetime');
    %opts=opts.setvaropts('date','InputFormat','yyyy_MM_dd');    
    opts = opts.setvaropts('acquisition_date','DatetimeFormat','dd-MMM-uuuu');    
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    histology_table=readtable(P.histology_table_path,opts);
    histology_table.probe_serial = [histology_table.probe_1 histology_table.probe_2];
    histology_table = removevars(histology_table,["probe_1","probe_2"]);
    histology_table.channels=arrayfun(@strsplit,histology_table.channels,'UniformOutput',false);
end