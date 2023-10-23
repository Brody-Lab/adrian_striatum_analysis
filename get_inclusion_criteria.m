function T = get_inclusion_criteria()
% RETURN a table containing the UberPhys unit inclusion criteria table
    P=get_parameters;
    path=P.unit_inclusion_criteria_table_path;
    if ~isfile(path)
        error('cannot find the inclusion criteria table')
    end
    opts = detectImportOptions(path,'TextType','string');
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    T=readtable(path,opts);
    warning('on','MATLAB:table:ModifiedAndSavedVarnames');
end