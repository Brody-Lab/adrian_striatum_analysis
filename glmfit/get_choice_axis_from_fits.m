function choice_axis = get_choice_axis_from_fits(stats,varargin)
   choice_axis = get_encoding_axis_from_fits(stats,{'cpoke_out_left','cpoke_out_right'},'time_lim_s',[-0.5 0],'time_average',true,'unit_norm',true,varargin{:}); 
end