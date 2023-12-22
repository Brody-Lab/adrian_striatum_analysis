function [fits_table,log] = select_cells_and_fits(glmfit_log,varargin)

    % creates a cells table with a fits column for all cells matching input
    % criteria and for the run matching the input criteria.
    % Takes all arguments that select_glmfit_runs and select_cells take.
    
    % this is the easiest way to just load up a bunch of fits
    
    % run `glmfit_log = load_glmfit_catalog()` first to load glmfit_log.
    % That's a big bottleneck.... Keep it in your workspace for later.

    P = get_parameters();
    log = select_glmfit_runs(glmfit_log,varargin{:}); % log is for only a single run
    if isempty(log)
        fits_table=[];
        return
    end
    cells_table = select_cells('is_in_dorsal_striatum',true,varargin{:});
    [fits,is_saved] = get_glm_fits(cells_table.recording_name,cells_table.cellno,log);
    fits_table = [cells_table(is_saved,:) fits];
    [~,idx] = sort(fits_table.AP);
    fits_table = fits_table(idx,:);
end