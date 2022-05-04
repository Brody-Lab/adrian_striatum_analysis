function fits = get_fits_for_session(recording_name,run,varargin)
    p=inputParser;
    p.addRequired('recording_name',@(x)validateattributes(x,{'string'},{'scalar'}));
    p.addParameter('cellno',[],@(x)validateattributes(x,{'numeric'},{'vector','positive','integer'}));
    p.addRequired('run',@(x)validateattributes(x,{'string'},{'scalar'}));
    p.parse(recording_name,run,varargin{:});
    params = p.Results;
    
    if isempty(params.cellno)
        cells_table = select_cells('recording_name',recording_name);            
        cellno = cells_table.cellno;
    end
    fits = get_glm_fits(recording_name,cellno,run,varargin{:});   
end