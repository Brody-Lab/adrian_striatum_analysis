function stats = add_click_diff_to_stats(stats)
    % adds pref click to stats struct as if it was a covariate
    %% extract relevant data from stats
    left_clicks = get_combined_weights(stats,'left_clicks');
    right_clicks = get_combined_weights(stats,'right_clicks');
    click_diff= right_clicks;
    click_diff.data = right_clicks.data - left_clicks.data;
    for i=numel(stats):-1:1
        stats(i).ws.click_diff = stats(i).ws.left_clicks;
        stats(i).ws.click_diff.data = click_diff.data(i,:);
    end
end