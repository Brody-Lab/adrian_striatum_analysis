function make_joined_cells_files()
    P=get_parameters();
    T=load_sessions_table;
    u=unique([T.rat,T.sessid],'rows');
    for i=1:size(u,1)
        x{i} = find(T.rat == u(i,1) & T.sessid == double(u(i,2)));
    end
    pairs=x(cellfun(@(x)numel(x)==2,x));
    joined_cells_path = fullfile(P.data_path,'cells','joined');
    if ~isfolder(joined_cells_path)
        mkdir(joined_cells_path);
    end
    for i=1:numel(pairs)
        fprintf('---\nMaking joined file for recordings %s and %s.\n',T.recording_name(pairs{i}(:)));
        Cells_1 = load_Cells_file(T.recording_name(pairs{i}(1)));
        Cells_2 = load_Cells_file(T.recording_name(pairs{i}(2)));
        Cells = join_cells_files(Cells_1,Cells_2);
        name = [char(T.rat(pairs{i}(1))),'_',datestr(T.sess_date(pairs{i}(1)),'yyyy_mm_dd')];
        if ~isfolder(fullfile(joined_cells_path,name))
            mkdir(fullfile(joined_cells_path,name));
        end
        filename = fullfile(joined_cells_path,name,[name,'_Cells.mat']);tic;
        fprintf('Saving joined file to %s...\n',filename);
        save(filename,'-struct','Cells');   
        fprintf(' took %s.\n',timestr(toc));
    end



end