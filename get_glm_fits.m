function [fits,is_responsive] = get_glm_fits(recording_name,cellno,run,varargin)

    % load fits given keys. together, recording_name and cellno uniquely
    % identify cells in the database. run specifies the fitting run and is
    % scalar. note that it takes about 0.1s to load each stats file over
    % the network so for many thousands of cells this can take a few
    % minutes.

    %% parse and validate inputs
    P=get_parameters;
    p=inputParser;
    p.addRequired('recording_name',@(x)validateattributes(x,{'string'},{'vector','nonempty'}));
    p.addRequired('cellno',@(x)validateattributes(x,{'numeric'},{'vector','nonempty','positive','integer'}));
    p.addRequired('run',@(x)validateattributes(x,{'string'},{'scalar','nonempty'}));
    p.parse(recording_name,cellno,run,varargin{:});    
    params = p.Results;
    n_cells = numel(recording_name);
    if numel(cellno)~=n_cells
        error('Number of elements of recording_name must match number of elements of cellno.');
    end
    
    %% add is_responsive column to cells_table based on run being used
    is_responsive = is_glmfit_responsive_cell(recording_name,cellno,run); % write this function
    recording_name = recording_name(is_responsive);
    cellno = cellno(is_responsive);
    n_cells = numel(recording_name);    
    
    %% load fits for each cell
    fits = table();
    for i=n_cells:-1:1
        fits.stats_path(i) = fullfile(P.fit_path,'fits',recording_name(i),run,'stats',['glmfit_stats_cell',num2str(cellno(i)),'.mat']);
        fits.stats(i) = load(fits.stats_path(i));
    end
    
end