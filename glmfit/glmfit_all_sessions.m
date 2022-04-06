function glmfit_all_sessions(varargin)
    % fits all sessions by assigning a cluster job to each session
    %% parse and validate inputs
    P=get_parameters;
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('kfold',1,@(x)validateattributes(x,{'numeric'},{'scalar','integer','>',0}));
    p.addParameter('bin_size_s',0.001,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));  % resolution of the model. predictions have this bin size.  
    p.addParameter('phi',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % value of C (adaptation state) at time lag 0
    p.addParameter('tau_phi',0,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % recovery time constant of C
    p.addParameter('fit_adaptation',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % whether or not to fit phi and tau_phi for each neuron using gradient descent
    p.addParameter('choice_time_back_s',0.75,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % choice kernels extend backwards acausally in time before stimulus end by this many seconds
    p.addParameter('time_per_job',5,@(x)validateattributes(x,{'numeric'},{'scalar','positive'})); % max time per job IN HOURS 
    p.addParameter('distribution','poisson',@(x)validateattributes(x,{'char'},{'nonempty'})); % max time per job IN HOURS     
    p.addParameter('include_mono_clicks',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('job_array',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % use a job array to parallelize over cells (useful if you are fitting adaptation params which takes a long time)
    p.parse(varargin{:});
    params=p.Results;    
    paths = get_data_paths('warn_existence',true,varargin{:});
    if P.on_tiger
        % if on cluster, run computationally expensive jobs first
        S = load_sessions_table();
        [~,order] = sort(S.n_clusters.*S.n_completed_acc_trials_no_stim,'descend');
        paths = paths(order);        
    end
    if ~all([paths.all_exist])
        error('Cannot continue because database is not fully accessible.');
    end
    time_string =datestr(now,'YYYY_mm_DD_HH_MM_SS');
    for i=1:length(paths)
        fprintf('\n---------Submitting job %g of %g-----------\n',i,length(paths));
        job_name = strcat(paths(i).recording_name,'_glmfit_',time_string);
        fprintf('Job name will be %s.\n',job_name);
        output_dir=fullfile(strrep(paths(i).parent_dir,'cells','fits'),['glmfit_',time_string]);
        mkdir(output_dir);
        fprintf('   Made output directory: %s\n   ',output_dir); 
        if params.job_array && P.on_tiger
            matlab_command = sprintf(['fit_glm_to_Cells(''%s'',''save_path'',''%s'',''save'',true,''bin_size_s'',%g,',...
                '''kfold'',%g,''fit_adaptation'',logical(%g),''phi'',%0.10g,''tau_phi'',%0.10g,''choice_time_back_s'',%0.10g,''distribution'',''%s'',''include_mono_clicks'',logical(%g));'],...
                paths(i).cells_file,output_dir,params.bin_size_s,params.kfold,params.fit_adaptation,params.phi,params.tau_phi,params.choice_time_back_s,params.distribution,params.include_mono_clicks);
            save_param_command = matlab_command(1:end-2);
            save_param_command=[save_param_command,[',''fit'',false);exit;']];
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
            sbatch_command = sprintf('sbatch -e %s -o %s -t %g -J "%s" --array=%s --mem-per-cpu=4G submit_matlab_job.slurm %s',error_file,out_file,round(params.time_per_job*60),job_name,array_string(1:end-1),matlab_command);
            system();   
        else
            error_file = fullfile(output_dir,'job%A.stderr');
            out_file = fullfile(output_dir,'job%A.stdout');             
            matlab_command = sprintf(['"fit_glm_to_Cells(''%s'',''save_path'',''%s'',''save'',true,''bin_size_s'',%g,',...
                '''kfold'',%g,''fit_adaptation'',logical(%g),''phi'',%0.10g,''tau_phi'',%0.10g,''choice_time_back_s'',%0.10g,''distribution'',''%s'',''include_mono_clicks'',logical(%g));exit;"'],...
                paths(i).cells_file,output_dir,params.bin_size_s,params.kfold,params.fit_adaptation,params.phi,params.tau_phi,params.choice_time_back_s,params.distribution,params.include_mono_clicks);            
            sbatch_command = sprintf('sbatch -e %s -o %s -t %g -J "%s" --mail-type=FAIL,TIME_LIMIT submit_matlab_job.slurm %s',error_file,out_file,round(params.time_per_job*60),job_name,matlab_command);
        end
        if P.on_tiger
            fprintf('Sending following system command to initiate job:\n   %s\n',sbatch_command);
            system(sbatch_command);                       
        else
            eval(matlab_command(2:end-1)); %remove double quotation marks in character sequence
        end
    end   
end