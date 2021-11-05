function T = add_missing_sessions_to_sessions_table()

% this is a hacky script to add sessions that are not yet sorted to the
% sessions table. it uses the fields of the recordings log and
% npx_penetration and leaves the result not filled in.

P  =get_parameters();
recordings_table = read_recordings_log(P.recordings_path);
T = load_sessions_table();
sz = height(T);
penetrations = npx_penetrations;

missing_sessions = find(recordings_table.striatum_glm==1 & ismissing(recordings_table.curated_cells_file) & ismissing(recordings_table.cells_file));

for i=missing_sessions(:)'
    T.sess_date(end+1) = recordings_table.date(i);
    T.rat{end} = char(recordings_table.rat_name(i));
    T.probe_serial(end) = recordings_table.probeSerialNumber(i);
    pen_idx = find(strcmp({penetrations.serial},num2str(T.probe_serial(end))) & strcmp({penetrations.rat},T.rat(end)));
    if length(pen_idx)~=1
        error('');
    end
    T.days_implanted(end) = days(T.sess_date(end) - penetrations(pen_idx).date_implanted);
    T.craniotomy_AP(end) = penetrations(pen_idx).craniotomy_AP;
    T.craniotomy_ML(end) = penetrations(pen_idx).craniotomy_ML;
    T.depth_inserted(end) = penetrations(pen_idx).depth_inserted;   
end

end



