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
    p.addParameter('fit_adaptation',true,@(x)validateattributes(x,{'logical'},{'scalar'})); % whether or not to fit phi and tau_phi for each neuron using gradient descent
    p.addParameter('choice_time_back_s',0.75,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % choice kernels extend backwards acausally in time before stimulus end by this many seconds
    p.addParameter('time_per_job',23.99,@(x)validateattributes(x,{'numeric'},{'scalar','positive'})); % max time per job IN HOURS 
    p.addParameter('include_mono_clicks',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('job_array',true,@(x)validateattributes(x,{'logical'},{'scalar'})); % use a job array to parallelize over cells (useful if you are fitting adaptation params which takes a long time)
    p.parse(varargin{:});
    params=p.Results;    
    cells_paths = get_data_paths('data_path',fullfile(P.data_path,'cells'),varargin{:});
    fits_paths = get_data_paths('data_path',fullfile(P.data_path,'fits'),'parent_dir',true,varargin{:},'must_exist',false);
    time_string =datestr(now,'YYYY_mm_DD_HH_MM_SS');
    for i=1:length(cells_paths)
        fprintf('\n---------Submitting job %g of %g-----------\n',i,length(cells_paths));
        if ~exist(cells_paths{i},'file') 
            error('Cells file not found: %s.',cells_paths{i});
        end
        fprintf('   Found data file at %s\n',cells_paths{i});
        try
            [rat,date] = extract_info_from_npx_path(cells_paths{i});
        catch
            warning('Failed to extract rat and date from path so job name will be uninformative.\n');
            rat='';
            date='';
        end
        output_dir=fullfile(fits_paths{i},['glmfit_',time_string]);
        mkdir(output_dir);
        fprintf('   Made output directory: %s\n   ',output_dir); 
        if params.job_array
            matlab_command = sprintf(['fit_glm_to_Cells(''%s'',''save_path'',''%s'',''save'',true,''bin_size_s'',%g,',...
                '''kfold'',%g,''fit_adaptation'',logical(%g),''phi'',%0.10g,''tau_phi'',%0.10g,''choice_time_back_s'',%0.10g,''include_mono_clicks'',logical(%g));'],...
                cells_paths{i},output_dir,params.bin_size_s,params.kfold,params.fit_adaptation,params.phi,params.tau_phi,params.choice_time_back_s,params.include_mono_clicks);
            save_param_command = matlab_command(1:end-2);
            save_param_command=[save_param_command,[',''fit'',false);']];
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
            matlab_command = ['"',matlab_command(1:end-2),[',''cellno'',id,''save_params'',false);'],'"'];
            system(sprintf('sbatch -e %s -o %s -t %g -J "%s" --array=%s submit_matlab_job.slurm %s',error_file,out_file,round(params.time_per_job*60),[rat,',',date,'_glm'],array_string(1:end-1),matlab_command));   
        else
            error_file = fullfile(output_dir,'job%A.stderr');
            out_file = fullfile(output_dir,'job%A.stdout');             
            matlab_command = sprintf(['"fit_glm_to_Cells(''%s'',''save_path'',''%s'',''save'',true,''bin_size_s'',%g,',...
                '''kfold'',%g,''fit_adaptation'',logical(%g),''phi'',%0.10g,''tau_phi'',%0.10g,''choice_time_back_s'',%0.10g,''include_mono_clicks'',logical(%g));"'],...
                cells_paths{i},output_dir,params.bin_size_s,params.kfold,params.fit_adaptation,params.phi,params.tau_phi,params.choice_time_back_s,params.include_mono_clicks);            
            system(sprintf('sbatch -e %s -o %s -t %g -J "%s" submit_matlab_job.slurm %s',error_file,out_file,round(params.time_per_job*60),[rat,',',date,'_glm'],matlab_command));               
        end
    end   
end