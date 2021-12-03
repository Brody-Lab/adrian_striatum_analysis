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
    p.parse(Cells,varargin{:});
    params = p.Results;
    
    % calculate stability and presence measures
    [stability,presence,mean_rate] = calculate_unit_stability(params.Cells);
    
    % select cells with higher stability, presence_ratio and meanrate than threshold
    units_select = is_in_region(Cells,params.region) & ...
                   stability < params.thresh_stability & ...
                   presence > params.thresh_presence & ...
                   mean_rate > params.thresh_meanrate;

end

function in_region = is_in_region(Cells,regions)
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