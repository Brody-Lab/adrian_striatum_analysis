function [scalar_proj,cos_sim,choice_norm,click_norm] = get_proj_click_choice(ws,varargin)
    choice_axis = get_choice_axis_from_fits(ws,varargin{:});  
    click_axis = get_click_diff_axis_from_fits(ws,varargin{:});
    [scalar_proj,cos_sim,choice_norm,click_norm] = get_projection_stats(choice_axis,click_axis,varargin{:});
end