function catalogue_glmfits(varargin)
    P=get_parameters;
    p=inputParser;
    p.addParameter('delete_empty',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.parse(varargin{:});
    params=p.Results;
    paths = get_cells_paths();
    for i=1:length(paths)
        parent_dir{i}=fileparts(paths{i});
        f{i} = dir(parent_dir{i});
        f{i} = {f{i}.name};
        f{i} = f{i}(contains(f{i},'glmfit_'));
        for k=1:length(f{i})
            stats_dir = fullfile(parent_dir{i},f{i}{k},'stats');
            if exist(stats_dir)
                tmp=dir(stats_dir);
                if params.delete_empty && sum([tmp.bytes])==0
                    % this run didn't saved out data
                    status=rmdir(fullfile(parent_dir{i},f{i}{k}),'s');
                    if ~status
                        error('Failed to delete empty folder: %s.\n',fullfile(parent_dir{i},f{i}{k}));
                    end
                else
                    params(k,i)=load(fullfile(parent_dir{i},f{i}{k},'params.mat'));
                end
            end
        end 
    end
    a=1;
end