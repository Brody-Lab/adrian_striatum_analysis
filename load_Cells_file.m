function Cells = load_Cells_file(recording_name,varargin)
    p=inputParser;
    p.addParameter('local',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.parse(varargin{:});
    params=p.Results;
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');    
    paths = get_data_paths('recording_name',recording_name,varargin{:});
    if params.local
        if isfile(paths.cells_file)
            fprintf('Loading Cells file from local copy @ %s ... \n',paths.cells_file);tic;
            Cells = load(paths.cells_file);
            fprintf(' took %s.\n-----------------\n',timestr(toc));      
        else
            fprintf('Local file does not exist yet. Loading Cells file from original location @ %s... \n',paths.original_cells_file);tic;
            Cells = load(paths.original_cells_file);
            fprintf('   ... took %s.\n-----------------\n',timestr(toc));   
        end
    else
        fprintf('Loading Cells file from original location @ %s... \n',paths.original_cells_file);tic;
        Cells = load(paths.original_cells_file);
        fprintf('   ... took %s.\n-----------------\n',timestr(toc));                 
    end 
    if isfield(Cells,'Cells')
        Cells=Cells.Cells;
    end
    warning('on','MATLAB:table:ModifiedAndSavedVarnames');    
end