function params = fix_responsive(params_path,cells_file_path)
    % loads a params file that may be subject to the bug where the
    % responsive cell list was not made properly, remakes the responsive
    % cell list, and saves back out the file.
    % cells_file_path can be a path or an actual cells structure
    if exist(params_path,'file')
        load(params_path);
    else
        error('File does not exist.');
    end
    if nargin<2
        cells_file_path = correct_file_path(params.dm.dspec.expt.param.mat_file_name);  
        cells_file_path = strrep(cells_file_path,'data','cells'); 
        if ~exist(cells_file_path,'file')
            error('Could not identify cells file: tried %s.',cells_file_path);
        end
    end
    fprintf('Fixing "responsive_enough" field in params file: %s.\n',params_path);
    [~,new_params] = fit_glm_to_Cells(cells_file_path,'save',false,'fit',false,'minResponsiveFrac',params.minResponsiveFrac,'minSpkParamRatio',params.minSpkParamRatio);
    params.responsive_enough = new_params.responsive_enough;
    save(params_path,'params');
end
