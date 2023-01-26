function agb_save_fig(varargin)
%agb_save_fig Save figure 
%
%  savefig(FILENAME) saves the current fig as .fig, .png, and .pdf in the directory
%  where the calling function lives (useful for saving output of day-by-day
%  analyses) with filename FILENAME (FILENAME must be a compact filename
%  without extension).
%
%  savefig(FILENAME,H) saves specified figure handle.
%
%  savefig(FILENAME, ... , PATH) saves to PATH instead of calling function
%  directory
%
%  savefig(...,'format',FORMAT) saves the file in the formats specified by
%  the cell array FORMAT (e.g. {'fig','pdf','bmp'})

% AGB 2023

st=dbstack('-completenames');

p=inputParser;
p.addRequired('filename',@(x)validateattributes(x,{'char','string'},{'nonempty'}));
p.addOptional('h',gcf,@(x)validateattributes(x,{'matlab.ui.Figure'},{'scalar','nonempty'}));
p.addOptional('path',fileparts(st(2).file),@(x)validateattributes(x,{'char','string'},{'nonempty'}));
p.addParameter('format',{'fig','pdf','png'},@(x)validateattributes(x,{'cell'},{'nonempty'}));
p.parse(varargin{:});
params=p.Results;

if ~isscalar(string(params.filename)) || isfile(params.filename)
    error('"FILENAME" must be convertible to a scalar string and should be in compact form (i.e. not full file path.');
end

if ~isscalar(string(params.path)) && isfolder(params.path)
    error('"PATH" must be convertible to a scalar string which is a valid folder.');
end

valid_formats = {'tif','tiff','pdf','jpeg','jpg','png','fig','eps'};
for i=1:numel(params.format)
    validatestring(params.format{i},valid_formats);
end

for i=1:numel(params.format)
    fullname = fullfile(params.path,[params.filename,'.',params.format{i}]);    
    fprintf('Saving %s.\n',fullname);
    if params.format{i}=="fig"
        savefig(params.h,fullname);
    else
        exportgraphics(params.h,fullname);
    end
end