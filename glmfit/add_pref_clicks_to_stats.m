function stats = add_pref_clicks_to_stats(stats,pref_mode)
    % adds pref click to stats struct as if it was a covariate
    %% extract relevant data from stats
    left_clicks = get_combined_weights(stats,'left_clicks');
    right_clicks = get_combined_weights(stats,'right_clicks');
    covariate_stats = {stats.covariate_stats};
    switch pref_mode 
        case "deviation"
            prefers_right = string( cellfun(@(x)x.pref_click_by_deviation,covariate_stats,'uni',0) )=="right";                           
            is_negative = cellfun(@(x)x.pref_click_max_deviation,covariate_stats)<1;            
        case "average" 
            prefers_right = string(cellfun(@(x)x.pref_click_by_average,covariate_stats,'uni',0))=="right";              
            is_negative = cellfun(@(x)x.pref_click_average,covariate_stats)<1;     
        case "both"
            prefers_right = string(cellfun(@(x)x.pref_click_by_deviation,covariate_stats,'uni',0))=="right";       
            prefers_right2 = string(cellfun(@(x)x.pref_click_by_average,covariate_stats,'uni',0))=="right";   
            neither = xor(prefers_right,prefers_right2);
            prefers_right = double(prefers_right);
            prefers_right(neither) = NaN;
            
            is_negative = cellfun(@(x)x.pref_click_average,covariate_stats)<1;  % arbitrary choice                  
    end  
    pref_clicks= right_clicks;
    pref_clicks.data(prefers_right==0,:) = left_clicks.data(prefers_right==0,:);
%    pref_clicks.data(is_negative==1,:) =   -pref_clicks.data(is_negative==1,:);    
    for i=numel(stats):-1:1
        stats(i).ws.pref_clicks = stats(i).ws.left_clicks;
        stats(i).ws.pref_clicks.data = pref_clicks.data(i,:);
    end
end