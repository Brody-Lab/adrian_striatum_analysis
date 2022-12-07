function units_select = select_units(Cells,varargin)
    % Pick which units to use based on region, stability, and presence
    % (based on Adrian definitions). Return boolean vector of units to take
    
    %   units_select = select_units(Cells): use default params to exclude units
    
    %   select_units(Cells,'thresh_stability', THRESH_VALUE): use a
    %   different stability value to exclude cells, default 5
    
    %   select_units(Cells,'thresh_presence', THRESH_VALUE): use a
    %   different presence on trials value to exclude cells, default 0.5  
    

    p=inputParser;
    p.KeepUnmatched=true;
    p.addRequired('Cells',@(x)validateattributes(x,{'struct'},{'scalar'}));
    p.addParameter('thresh_stability',5,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); 
    p.addParameter('thresh_presence',0.5,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('thresh_meanrate',0.1,@(x)validateattributes(x,{'numeric'},{'nonnegative','scalar'}));    
    p.addParameter('region',0:100,@(x)validateattributes(x,{'numeric','char','cell'},{'nonempty'}));
    ncells=numel(Cells.raw_spike_time_s);
    p.addParameter('exclude_cells',false(ncells,1),@(x)validateattributes(x,{'logical'},{'vector'}));     
    p.parse(Cells,varargin{:});
    params = p.Results;
    
    % calculate stability and presence measures
    [stability,presence,mean_rate] = calculate_unit_stability(params.Cells);
    
    % select cells with higher stability, presence_ratio and meanrate than threshold
    stable = stability < params.thresh_stability;
    present = presence > params.thresh_presence;
    good_rate = mean_rate > params.thresh_meanrate;
    
    units_select = stable & present & good_rate & ~params.exclude_cells(:);   
    
    fprintf('-----\nUser excluded %g of %g (%d%%) cells.\n-----\n',sum(params.exclude_cells),ncells,round(mean(params.exclude_cells)*100));
    
    fprintf(['%g of %g (%d%%) cells passed stability threshold (%g).\n',...
        '%g of %g (%d%%) cells passed presence threshold (trial fraction = %g).\n',...
        '%g of %g (%d%%) cells passed rate threshold (%g sp/s).\n------\n',...
        '%g of %g cells (%d%%) selected.'],...
        sum(stable),ncells,round(100*mean(stable)),params.thresh_stability,sum(present),ncells,round(100*mean(present)),params.thresh_presence,...
        sum(good_rate),ncells,round(100*mean(good_rate)),params.thresh_meanrate,sum(units_select),ncells,round(100*mean(units_select)));
               
    if ~ismember('region',p.UsingDefaults)
        units_select = is_in_region(Cells,params.region) & units_select;
    end
    
end

function in_region = is_in_region(Cells,regions)
    Cells.regions=Cells.regions(:);
    if isnumeric(regions)
        in_region = ismember(Cells.regions,regions);
        return
    elseif ischar(regions)
        for i=1:length(Cells.penetration.regions)
            match(i) = any(strcmp(regions,Cells.penetration.regions(1).name));
        end
        if sum(match)==1
            in_region = is_in_region(Cells,find(match));
            return
        else
            error('None or multiple region matches in Cells.penetration.regions');
        end
    elseif iscell(regions)
        in_region = false(size(Cells.regions));
        for i=1:length(regions)
            in_region = in_region | is_in_region(Cells,regions{i});
        end
    end
end