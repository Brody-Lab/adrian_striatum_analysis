function filenames = plot_glm_covariates(recording_name,cellno,log)
    n_cells = numel(cellno);
    if isscalar(recording_name)
        recording_name = repmat(recording_name,n_cells,1);
    end
    [fits,is_saved] = get_glm_fits(recording_name,cellno,log);
    not_saved = cellno(~is_saved);
    saved_cells = cellno(is_saved);
    not_saved_recording_name = recording_name(~is_saved);
    for i=1:numel(not_saved)
        fprintf('No fit found for %s cell %d.\n',not_saved_recording_name(i),not_saved(i));    
    end
    parfor i=1:numel(saved_cells)
        h(i)=figure;           
        plotGLM.plotFittedCovariates(fits.stats(i));
        filenames(i)=strrep(fits.stats_path(i),'stats','plot_covariates');
        filenames(i)=strrep(filenames(i),'.mat','.png');            
        [plot_fldr,~] = fileparts(filenames(i));
        if ~isfolder(plot_fldr)
            mkdir(plot_fldr);
        end
        saveas(h(i),filenames(i));
        close(h(i));
    end
end