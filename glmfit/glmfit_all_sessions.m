function glmfit_all_sessions(varargin)
    % fits all sessions by assigning a cluster job to each session
    %% parse and validate inputs
    P=get_parameters;
    paths = get_data_paths('warn_existence',true,varargin{:});
    
    % params exclusive to this function
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('job_array',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % use a job array to parallelize over cells (useful if you are fitting adaptation params which takes a long time)
    p.addParameter('time_per_job',4,@(x)validateattributes(x,{'numeric'},{'scalar','positive'})); % max time per job IN HOURS     
    p.addParameter('partition','all',@(x)validateattributes(x,{'char','string'},{'nonempty'}));
    p.addParameter('requeue',false,@(x)validateattributes(x,{'logical'},{'scalar'})); 
    p.addParameter('alpha',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % elastic net parameter (0=ridge, 1=lasso) 
    p.addParameter('glmnet_thresh',1e-6,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));    
    p.addParameter('sbatch_retry_frequency_mins',5,@(x)validateattributes(x,{'numeric'},{'scalar','positive'})); % retry sbatch submission after this many minutes
    p.addParameter('jobid',datestr(now,'YYYY_mm_DD_HH_MM_SS'),@(x)validateattributes(x,{'char','string'},{'nonempty'}));
    p.addParameter('paths',1:numel(paths),@(x)validateattributes(x,{'numeric'},{'integer','increasing','nonempty','positive','<=',numel(paths)})); %indices of sessions to run (as they would run after sorting step)
    p.parse(varargin{:});    
    
    %glm fitting params
    params = glmfit_parse_params(varargin{:});
    
    try
        datetime(p.Results.jobid,'InputFormat','yyyy_MM_dd_HH_mm_ss');
    catch
        error('job id %s could not be parsed as a date string.',p.Results.jobid);
    end
    
    if P.on_cluster
        % if on cluster, run computationally expensive jobs first
        S = load_sessions_table();
        if p.Results.job_array
            [~,order] = sort(S.n_trials,'descend');            
        else
            [~,order] = sort(S.n_clusters.*S.n_trials,'descend');
        end
        paths = paths(order);        
    end
    if ~all([paths.all_exist])
        error('Cannot continue because database is not fully accessible.');
    end
    for i=p.Results.paths(:)'
        fprintf('\n---------Submitting job %g of %g-----------\n',i,length(paths));
        job_name = strcat(paths(i).recording_name,'_glmfit_',p.Results.jobid);
        fprintf('Job name will be %s.\n',job_name);
        output_dir=fullfile(strrep(paths(i).parent_dir,'cells','fits'),['glmfit_',p.Results.jobid]);
        if isfolder(output_dir)
            fprintf('   Output directory already exists: %s\n   ',output_dir);             
        else
            mkdir(output_dir);
            fprintf('   Made output directory: %s\n   ',output_dir); 
        end
        if p.Results.job_array && P.on_cluster
            matlab_command = sprintf(['fit_glm_to_Cells(''%s'',''save_path'',''%s'',''save'',true,''bin_size_s'',%0.10g,',...
                '''kfold'',%g,''phi'',%0.10g,''tau_phi'',%0.10g,''choice_time_back_s'',%0.10g,''alpha'',%0.10g,''glmnet_thresh'',%0.10g,',...
                '''distribution'',''%s'',''include_mono_clicks'',logical(%g),''use_trial_history'',logical(%g),''nClickBins'',%g,',...
                '''separate_clicks_by'',''%s'',''separate_clicks_by_outcome'',logical(%g));'],...
                paths(i).cells_file,output_dir,params.bin_size_s,params.kfold,params.phi,params.tau_phi,params.choice_time_back_s,...
                params.alpha,params.glmnet_thresh,params.distribution,params.include_mono_clicks,params.use_trial_history,...
                params.nClickBins,params.separate_clicks_by,params.separate_clicks_by_outcome);
            save_param_command=[matlab_command(1:end-2),',''fit'',false);'];
            fprintf('Running fit_glm_to_Cells to save params before batch fitting.\n');
            eval(save_param_command);
            params_path=fullfile(output_dir,'glmfit_params.mat');
            if exist(params_path,'file')
               glmfit_params=load(params_path);
            else
               error('Could not load saved params file or it was not saved: %s\n',params_path); 
            end
            error_file = fullfile(output_dir,'job%A_cell%a.stderr');
            out_file = fullfile(output_dir,'job%A_cell%a.stdout');            
            array_string=sprintf('%g,',glmfit_params.params.cellno(glmfit_params.params.responsive_enough));
            matlab_command = ['"',matlab_command(1:end-2),',''cellno'',id,''save_params'',false);','exit;"'];
            if p.Results.requeue
                sbatch_command = sprintf('sbatch -e %s -o %s -t %g -J "%s" --array=%s --mem-per-cpu=5G -p %s --requeue submit_matlab_job.slurm %s',...
                    error_file,out_file,round(p.Results.time_per_job*60),job_name,array_string(1:end-1),p.Results.partition,matlab_command);
            else
                sbatch_command = sprintf('sbatch -e %s -o %s -t %g -J "%s" --array=%s --mem-per-cpu=5G -p %s submit_matlab_job.slurm %s',...
                    error_file,out_file,round(p.Results.time_per_job*60),job_name,array_string(1:end-1),p.Results.partition,matlab_command);                
            end
        else
            error_file = fullfile(output_dir,'job%A.stderr');
            out_file = fullfile(output_dir,'job%A.stdout');             
            matlab_command = sprintf(['"fit_glm_to_Cells(''%s'',''save_path'',''%s'',''save'',true,''bin_size_s'',%g,',...
                '''kfold'',%g,''phi'',%0.10g,''tau_phi'',%0.10g,''choice_time_back_s'',%0.10g,''alpha'',%0.10g,''glmnet_thresh'',%0.10g,'...
                '''distribution'',''%s'',''include_mono_clicks'',logical(%g),''use_trial_history'',logical(%g),''nClickBins'',%g,',...
                '''separate_clicks_by'',''%s'',''separate_clicks_by_outcome'',logical(%g));exit;"'],...
                paths(i).cells_file,output_dir,params.bin_size_s,params.kfold,params.phi,params.tau_phi,params.choice_time_back_s,params.alpha,...
                params.glmnet_thresh,params.distribution,params.include_mono_clicks,params.use_trial_history,params.nClickBins,...
                params.separate_clicks_by,params.separate_clicks_by_outcome);       
            if p.Results.requeue
                sbatch_command = sprintf('sbatch -e %s -o %s -t %g -J "%s" --mail-type=FAIL,TIME_LIMIT -p %s --requeue submit_matlab_job.slurm %s',...
                    error_file,out_file,round(p.Results.time_per_job*60),job_name,p.Results.partition,matlab_command);
            else
                sbatch_command = sprintf('sbatch -e %s -o %s -t %g -J "%s" --mail-type=FAIL,TIME_LIMIT -p %s submit_matlab_job.slurm %s',...
                    error_file,out_file,round(p.Results.time_per_job*60),job_name,p.Results.partition,matlab_command);                
            end
        end
        if P.on_cluster
            fprintf('Sending following system command to initiate job:\n   %s\n',sbatch_command);
            while true
                output = system(sbatch_command);    
                if output
                    fprintf('sbatch command failed. Trying again in %d minutes.\n',round(p.Results.sbatch_retry_frequency_mins));
                    pause(p.Results.sbatch_retry_frequency_mins*60);
                else
                    break
                end
            end
                
        else
            eval(matlab_command(2:end-1)); %remove double quotation marks in character sequence
        end
    end   
end