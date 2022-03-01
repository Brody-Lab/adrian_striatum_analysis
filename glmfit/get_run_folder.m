function run_folders = get_run_folder(runs,recording_names)
    P=get_parameters;
    p=inputParser;
    p.addRequired('run',@(x)validateattributes(x,{'string'},{'vector','nonempty'}));    
    p.addRequired('recording_names',@(x)validateattributes(x,{'string'},{'vector','nonempty'}));
    p.parse(runs,recording_names);    
    n=numel(runs);
    if n~=numel(recording_names)
        if n==1
            runs = repmat(runs,numel(recording_names),1);
            n=numel(recording_names);
        else
            error('runs and recording_names must be vectors of same length or runs must be scalar.'); 
        end
    end
    for i=n:-1:1
       run_folders(i,1) = string(fullfile(P.fit_path,'fits',recording_names(i),runs(i)));
    end
end
