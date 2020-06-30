function paths=get_tigress_cells_paths(varargin)
    P=get_parameters;
    p=inputParser;
    p.addParameter('rat','',@(x)validateattributes(x,{'char','cell'},{}));
    p.parse(varargin{:});
    params=p.Results;
    recordings_table = read_recordings_log(P.recordings_path);
    curated_cells_files = recordings_table.curated_cells_file(recordings_table.striatum_glm==1);
    cells_files = recordings_table.cells_file(recordings_table.striatum_glm==1);
    rat_names = recordings_table.rat_name(recordings_table.striatum_glm==1);
    cells_files(~ismissing(curated_cells_files))=curated_cells_files(~ismissing(curated_cells_files));
    fix_path = @(x)strrep(strrep(strrep(char(x),'"',''),'\',filesep),'/',filesep);
    for i=1:length(cells_files)
        paths{i}=fullfile(P.tiger_volume,'data',fix_path(char(regexprep(cells_files(i),'.*Adrian(.*)','$1'))));       
        if ~exist(paths{i},'file')
            %error('File not found: %s.',paths{i});
        end
    end
    if ~isempty(params.rat)
        if ischar(params.rat)
            params.rat={params.rat};
        end
        paths = paths(ismember(rat_names,params.rat));
    end
end