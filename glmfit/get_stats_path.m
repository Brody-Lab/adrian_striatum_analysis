function stats_path = get_stats_path(runs,recording_names,cellno)
    p=inputParser;
    p.addRequired('run',@(x)validateattributes(x,{'string'},{'vector','nonempty'}));    
    p.addRequired('recording_names',@(x)validateattributes(x,{'string'},{'vector','nonempty'}));
    p.addRequired('cellno',@(x)validateattributes(x,{'numeric'},{'vector','nonempty','positive'})); 
    p.parse(runs,recording_names,cellno);       
    if any( size(recording_names) ~= size(cellno))
        error('Recording names and cellnos must be vectors of the same size');
    end
   run_folders = get_run_folder(runs,recording_names);
   for i=numel(cellno):-1:1
       stats_path(i,1) = fullfile(run_folders(i),'stats',['glmfit_stats_cell',num2str(cellno(i)),'.mat']);
   end
end