function create_database(varargin)
    % moves files from bucket into the repo database, and creates index
    % files for easy searching of what's there
    P=get_parameters;
    p=inputParser;
    p.addParameter('make_cell_info',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % make cell_info and session_info files?
    p.addParameter('remake',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % add only new sessions?
    p.addParameter('recopy',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % write over cells files?   
    p.addParameter('resave',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % resave cells files after adding fields?
    p.parse(varargin{:});
    params = p.Results;
    recordings_table = read_recordings_log(P.recordings_path);
    tagged = recordings_table.D2Phototagging==1;
    tagged = tagged(recordings_table.striatum_glm==1);
    laser_power_mW = recordings_table.laser_power_mW;
    laser_power_mW = laser_power_mW(recordings_table.striatum_glm==1);    
    session_notes = recordings_table.notes;
    session_notes = session_notes(recordings_table.striatum_glm==1);    
    D2Phototagging = recordings_table.D2Phototagging;
    D2Phototagging = D2Phototagging(recordings_table.striatum_glm==1);        
    curated_cells_files = recordings_table.curated_cells_file(recordings_table.striatum_glm==1);
    cells_files = recordings_table.cells_file(recordings_table.striatum_glm==1);
    cells_files(~ismissing(curated_cells_files))=curated_cells_files(~ismissing(curated_cells_files));
    if any(ismissing(cells_files))
        warning('%g missing cells files. Skipping.\n',sum(ismissing(cells_files)));
    end    
    tagged = tagged(~ismissing(cells_files));
    laser_power_mW = laser_power_mW(~ismissing(cells_files));
    session_notes = session_notes(~ismissing(cells_files));
    D2Phototagging = D2Phototagging(~ismissing(cells_files));    
    cells_files = cells_files(~ismissing(cells_files));
    fix_path = @(x)strrep(char(x),'"','');
    for i=1:length(cells_files)
%         if ~tagged(i)
%             continue
%         end
        destination=fullfile(P.data_path,'cells',char(regexprep(cells_files(i),'.*Adrian(.*)','$1')));       
        if ~exist(fix_path(cells_files(i)),'file')
            error('File not found: %s.',fix_path(cells_files(i)));
        end
        if ~isdir(fileparts(destination))
            mkdir(fileparts(destination));
        end
        if ~isfile(fix_path(destination)) || params.recopy 
            fprintf('Copying %s ----->\n   to %s ...',fix_path(cells_files(i)),fix_path(destination));tic
            copyfile(fix_path(cells_files(i)),fix_path(destination));
            fprintf(' took %s.\n',timestr(toc));   
        else
            fprintf('%s exists.\n',fix_path(destination));          
        end
        if params.make_cell_info
            [parent,~,~] = fileparts(fix_path(destination));
            cell_info_path = fullfile(parent,'cell_info.mat');    
            session_info_path = fullfile(parent,'session_info.mat');
            if isfile(cell_info_path) && ~params.remake 
                fprintf('%s exists.\n-----------------\n',cell_info_path);                         
            else
                fprintf('Loading Cells file  ...');tic;
                Cells = load(fix_path(destination));
                fprintf(' took %s.\n-----------------\n',timestr(toc));                
                fields=fieldnames(Cells);
                if length(fields)==1
                    Cells=Cells.(fields{1}); % if not saved with -struct flag
                end   
                Cells.laser_power_mW = laser_power_mW(i);
                Cells.session_notes = session_notes(i);
                Cells.D2Phototagging = D2Phototagging(i);
                Cells.mat_file_name = string(cells_files{i});
                if ~isfield(Cells,'ap_meta') && isfield(Cells,'meta')
                    Cells.ap_meta = Cells.meta.ap_meta;
                end    
                Cells = import_penetration(Cells);    
                Cells.n_clusters = numel(Cells.spike_time_s.cpoke_in);                
                if ~isfield(Cells,'ks_good')
                    % calculate refractory period violations a la Kilosort
                    for k=1:Cells.n_clusters
                        Cells.ks_good(k) = is_ks_good(Cells.raw_spike_time_s{k});
                    end
                end    
                try
                    Cells.probe_serial = Cells.ap_meta.imDatPrb_sn;
                catch
                    Cells.probe_serial = Cells.rec.ap_meta.imDatPrb_sn;        
                end
                time_to_clicks = Cells.Trials.stateTimes.clicks_on - Cells.Trials.stateTimes.cpoke_in;
                exclude_trials = Cells.Trials.violated | time_to_clicks<0.5 | Cells.Trials.laser.isOn;
                [Cells.autocorr,~,Cells.autocorr_fr_hz] = timescales.get_autocorr_from_Cells(Cells,'cpoke_in',[0 0.5],'exclude_trials',exclude_trials);                  
                if isfield(Cells,'sess_date')
                    try
                        Cells.sess_date = datetime(Cells.sess_date);
                    catch
                        Cells.sess_date = datetime(Cells.sess_date,'InputFormat','uuuu_MM_dd');
                    end
                elseif unique(Cells.sessid)==701531 % special case for A242 session
                    Cells.sess_date = datetime('2019-06-03');
                elseif unique(Cells.sessid)==702016 % special case for A242 session
                    Cells.sess_date = datetime('2019-06-06');      
                elseif unique(Cells.sessid)==703121 % special case for A242 session     
                    Cells.sess_date = datetime('2019-06-10');              
                end
                Cells.rat = Cells.rat(1,:);                                
                Cells.sess_date = Cells.sess_date(1,:);                                
                Cells.days_implanted = days(datetime(Cells.sess_date) - datetime(Cells.penetration.date_implanted)); 
                if ~isfield(Cells,'last_modified') && isfield(Cells,'meta')
                    Cells.last_modified = Cells.meta.last_modified;
                end
                Cells.last_modified = datetime(Cells.last_modified);
                Cells.probe_serial = string(Cells.probe_serial);     
                if isfield(Cells,'PETH')
                    Cells = rmfield(Cells,'PETH');
                end         
                fprintf('Making cell_info  ...');tic;      
                if tagged(i)
                    Cells = tagging.compute_laser_modulation(Cells);
                end    
                Cells.sessid = Cells.sessid(1);
                Cells.hemisphere = string(Cells.penetration.hemisphere);
                cell_info = make_cell_info(Cells,tagged(i));
                fprintf(' took %s.\n-----------------\n',timestr(toc));                
                if params.resave
                    fprintf('saving Cells file  ...');tic;                    
                    save(fix_path(destination),'-struct','Cells');
                    fprintf(' took %s.\n-----------------\n',timestr(toc));                                    
                end
                save(cell_info_path,'cell_info');
                % make and save session info too
                session_info_fields = {'sess_date','rat','sessid','mat_file_name','probe_serial','days_implanted','n_clusters','craniotomy_AP','craniotomy_ML','depth_inserted','hemisphere','region_names','n_trials','violation_rate','percent_correct','laser_power_mW','D2Phototagging','session_notes'};
                for f=1:length(session_info_fields)
                    if isfield(Cells,session_info_fields{f})
                        session_info.(session_info_fields{f}) = Cells.(session_info_fields{f}); 
                    else
                        if session_info_fields{f}=="region_names"
                            session_info.region_names  = {Cells.penetration.regions.name};
                        elseif session_info_fields{f}=="n_trials"
                            session_info.n_trials = numel(Cells.Trials.is_hit);
                        elseif session_info_fields{f}=="violation_rate"
                            session_info.violation_rate = mean(Cells.Trials.violated);                
                        elseif session_info_fields{f}=="percent_correct"
                            session_info.percent_correct = 100*mean(Cells.Trials.is_hit(~Cells.Trials.violated));
                        else
                            session_info.(session_info_fields{f}) = Cells.penetration.(session_info_fields{f}); 
                        end
                    end
                end
                save(session_info_path,'session_info');                      
            end      
        end
    end
end