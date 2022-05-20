function [scalar_proj,scalar_proj_orth,cos_sim,choice_norm,click_norm] = get_proj_click_choice(ws,varargin)
    if nargin>1
        ws_choice = varargin{1};
        varargin = varargin(2:end);
    else
        ws_choice = ws;
    end
    choice_axis = get_choice_axis_from_fits(ws_choice,varargin{:});  
    click_axis = get_click_diff_axis_from_fits(ws,varargin{:});
    [scalar_proj,scalar_proj_orth,cos_sim,choice_norm,click_norm] = get_projection_stats(choice_axis,click_axis,varargin{:});
end