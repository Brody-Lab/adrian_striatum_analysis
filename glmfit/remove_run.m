function remove_run(run_folder)
    fprintf('Removing run folder %s  ...',run_folder);tic;
    if isfolder(run_folder)
        status=rmdir(run_folder,'s');
        if ~status
            error('Failed to delete empty folder: %s.\n',run_folder);
        end    
    else
        error('%s is not a valid folder.\n',run_folder);        
    end
    fprintf('took %s.\n',timestr(toc));
end