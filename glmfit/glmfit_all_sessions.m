function glmfit_all_sessions(varargin)
    % fits all sessions by assigning a cluster job to each session
    %% parse and validate inputs
    P=get_parameters;
    paths = get_data_paths('warn_existence',true,varargin{:});
    
    % params exclusive to this function
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('job_array',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % use a job array to parallelize over cells (useful if you are fitting adaptation params which takes a long time)
    p.addParameter('time_per_job',10,@(x)validateattributes(x,{'numeric'},{'scalar','positive'})); % max time per job IN HOURS     
    p.addParameter('partition','',@(x)validateattributes(x,{'char','string'},{}));
    p.addParameter('sbatch_retry_frequency_mins',5,@(x)validateattributes(x,{'numeric'},{'scalar','positive'})); % retry sbatch submission after this many minutes
    p.addParameter('jobid',datestr(now,'YYYY_mm_DD_HH_MM_SS'),@(x)validateattributes(x,{'char','string'},{'nonempty'}));
    p.addParameter('paths',1:numel(paths),@(x)validateattributes(x,{'numeric'},{'integer','increasing','nonempty','positive','<=',numel(paths)})); %indices of sessions to run (as they would run after sorting step)
    p.addParameter('ncpus',1,@(x)validateattributes(x,{'numeric'},{'scalar','positive','integer'}));    
    p.addParameter('ngpus',0,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','integer'}));  
    p.addParameter('memory_gb',7,@(x)validateattributes(x,{'numeric'},{'scalar','positive','integer'}));        
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
        C = load_cells_table();
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
        matlab_command = sprintf(['"fit_glm_to_Cells(''%s'',''save_path'',''%s'',''save'',true,''bin_size_s'',%0.10g,',...
            '''kfold'',%g,''phi'',%0.10g,''tau_phi'',%0.10g,''choice_time_back_s'',%0.10g,''alpha'',%0.10g,''glmnet_thresh'',%0.10g,',...
            '''distribution'',''%s'',''include_mono_clicks'',logical(%g),''use_trial_history'',logical(%g),''nClickBins'',%g,',...
            '''separate_clicks_by'',''%s'',''separate_clicks_by_outcome'',logical(%g),''fit_method'',''%s'',''parallel'',logical(%g),',...
            '''useGPU'',logical(%g),''relaxed'',logical(%g),''lambda'',%0.10g,''clicks_boundary'',{''%s'',''%s''},''separate_cpoke_out_by_choice'',logical(%g),',...
            '''clicks_boundary_offset'',[%0.10g %0.10g],''choice'',''%s'');"'],...
            paths(i).cells_file,output_dir,params.bin_size_s,params.kfold,params.phi,params.tau_phi,params.choice_time_back_s,...
            params.alpha,params.glmnet_thresh,params.distribution,params.include_mono_clicks,params.use_trial_history,...
            params.nClickBins,params.separate_clicks_by,params.separate_clicks_by_outcome,params.fit_method,params.parallel,params.useGPU,...
            params.relaxed,params.lambda,params.clicks_boundary{1},params.clicks_boundary{2},params.separate_cpoke_out_by_choice,...
            params.clicks_boundary_offset(1),params.clicks_boundary_offset(2),params.choice);        
        if p.Results.job_array && P.on_cluster
            save_param_command=[matlab_command(2:end-3),',''fit'',false);'];
            fprintf('Running fit_glm_to_Cells to save params before batch fitting.\n');
            eval(save_param_command);
            params_path=fullfile(output_dir,'glmfit_params.mat');
            if exist(params_path,'file')
               glmfit_params=load(params_path);
            else
               error('Could not load saved params file or it was not saved: %s\n',params_path); 
            end
            thisC=C(C.recording_name==paths(i).recording_name,:);        
            cellnos = find(apply_inclusion_criteria(thisC));            
            error_file = fullfile(output_dir,'job%A_cell%a.stderr');
            out_file = fullfile(output_dir,'job%A_cell%a.stdout');            
            array_string=sprintf('%g,',intersect(cellnos,glmfit_params.params.cellno(glmfit_params.params.responsive_enough)));
            matlab_command = ['"',matlab_command(2:end-3),',''cellno'',id,''save_params'',false);','exit;"'];                
        else
            error_file = fullfile(output_dir,'job%A.stderr');
            out_file = fullfile(output_dir,'job%A.stdout'); 
            array_string='0,';         
        end
        if P.on_cluster
            sbatch_command = sprintf('sbatch -e %s -o %s -t %g -J "%s" --array=%s -p %s --cpus-per-task=%d --gres=gpu:%d --mem-per-cpu=%dG submit_matlab_job.slurm %s',...
                error_file,out_file,round(p.Results.time_per_job*60),job_name,array_string(1:end-1),p.Results.partition,p.Results.ncpus,p.Results.ngpus,p.Results.memory_gb,matlab_command);            
            sbatch_command = strrep(sbatch_command,'--gres=gpu:0 ','');            
            sbatch_command = strrep(sbatch_command,'-p  ','');                        
            fprintf('Sending following system command to initiate job:\n   %s\n',sbatch_command);
            while true
                output=0;
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