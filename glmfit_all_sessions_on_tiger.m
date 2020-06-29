function glmfit_all_sessions_on_tiger(varargin)
    %% parse and validate inputs
    P=get_parameters;
    p=inputParser;
    p.addParameter('kfold',1,@(x)validateattributes(x,{'numeric'},{'scalar','integer','>',0}));
    p.addParameter('bin_size_s',0.001,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));  % resolution of the model. predictions have this bin size.  
    p.addParameter('phi',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % value of C (adaptation state) at time lag 0
    p.addParameter('tau_phi',0,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % recovery time constant of C
    p.addParameter('fit_adaptation',true,@(x)validateattributes(x,{'logical'},{'scalar'})); % whether or not to fit phi and tau_phi for each neuron using gradient descent
    p.addParameter('choice_time_back_s',0.75,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % choice kernels extend backwards acausally in time before stimulus end by this many seconds
    p.addParameter('time_per_job',23.99,@(x)validateattributes(x,{'numeric'},{'scalar','positive'})); % max time per job IN HOURS 
    p.addParameter('include_mono_clicks',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.parse(varargin{:});
    params=p.Results;    
    paths = get_tigress_cells_paths();
    time_string =datestr(now,'YYYY_mm_DD_HH_MM_SS');
    for i=1:length(paths)
        output_dir=fullfile(fileparts(paths{i}),['glmfit_',time_string]);
        error_file = fullfile(output_dir,'slurm.stderr');
        out_file = fullfile(output_dir,'slurm.stdout');       
        matlab_command = sprintf(['"fit_glm_to_Cells(''%s'',''output_dir'',''%s'',''save'',true,''bin_size_s'',%g,',...
            '''kfold'',%g,''fit_adaptation'',logical(%g),''phi'',%0.10g,''tau_phi'',%0.10g,''choice_time_back_s'',%0.10g,''include_mono_clicks'',logical(%g));"'],...
            paths{i},output_dir,params.bin_size_s,params.kfold,params.fit_adaptation,params.phi,params.tau_phi,params.choice_time_back_s,params.include_mono_clicks);
        system(sprintf('sbatch -e %s -o %s -t %g -J "striatum_glm" submit_matlab_job.slurm %s',error_file,out_file,round(params.time_per_job*60),matlab_command));   
    end   
end


