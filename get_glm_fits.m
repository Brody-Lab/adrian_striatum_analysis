function [fits,is_responsive] = get_glm_fits(recording_name,cellno,run,varargin)

    % load fits given keys. together, recording_name and cellno uniquely
    % identify cells in the database. run specifies the fitting run and is
    % scalar. note that it takes about 0.1s to load each stats file over
    % the network so for many thousands of cells this can take a few
    % minutes.

    %% parse and validate inputs
    p=inputParser;
    p.addRequired('recording_name',@(x)validateattributes(x,{'string'},{'vector','nonempty'}));
    p.addRequired('cellno',@(x)validateattributes(x,{'numeric'},{'vector','nonempty','positive','integer'}));
    p.addRequired('run',@(x)validateattributes(x,{'string','table'},{'nonempty'}));
    p.addParameter('parallel_load',true);
    p.parse(recording_name,cellno,run,varargin{:});
    params = p.Results;
    n_cells = numel(cellno);
    if istable(run) && ismember('run',run.Properties.VariableNames) && numel(unique(run.run))==1
        run=run.run(1);
    elseif isscalar(run) && isstring(run)
    else
        error('Run must be a glmfit_log table with only entries having the same run name or it must be a scalar string specifying the runname.');
    end        
    if numel(recording_name)~=n_cells && ~isscalar(recording_name)
        error('Number of elements of recording_name must match number of elements of cellno.');
    end
    if isscalar(recording_name)
        recording_name = repmat(recording_name,n_cells,1);
    end
    
    %% add is_responsive column to cells_table based on run being used
    is_responsive = is_glmfit_responsive_cell(recording_name,cellno,run);
    recording_name = recording_name(is_responsive);
    cellno = cellno(is_responsive);
    n_cells = numel(recording_name);    
    
    %% load fits for each cell
    fits = table();
    fprintf('get_glm_fits: Getting paths for stats files ');        tic;
    fits.stats_path = get_stats_path(run,recording_name,cellno);
    fprintf(' ... took %s.\n',timestr(toc));
    fprintf('get_glm_fits: Loading %g stats files ... ',n_cells);tic;  
    if params.parallel_load
        stats_path = fits.stats_path;        
        p = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(p)
            parpool(6);
        end    
        parfor i=1:n_cells % with 6 workers, parallel loading is about 4x faster than serial loading for fits on jukebox.
            stats(i,1) = load(stats_path(i));
        end
        fits.stats=stats;        
    else
        count=0;
        for i=n_cells:-1:1 
            count=count+1;
            fits.stats(i) = load(fits.stats_path(i));
            if mod(count,100)==0
                fprintf('Finished %g loads in %s.\n',count,timestr(toc));tic;
            end
        end
    end
    fprintf(' took %s.\n',timestr(toc));
end