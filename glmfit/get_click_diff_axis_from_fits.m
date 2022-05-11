function choice_axis = get_click_diff_axis_from_fits(ws,varargin)
   choice_axis = get_encoding_axis_from_fits(ws,{'left_clicks','right_clicks'},'time_lim_s',[0 1],'time_average',false,'unit_norm',false,varargin{:}); 
end