% GET_PARAMETERS Returns a structure with the parameters for analyses, such as file paths
%
%=OUTPUT
%
%   P
%       A structure with parameters for analyses
%
function P = get_parameters()
%% paths
P.repository_path = fileparts(mfilename('fullpath'));
P.recordings_path = fullfile(P.repository_path,'recordings_log.csv');
P.tiger_volume = '/tigress';

%% plotting
P.figure_image_format = {'png'};
P.figure_position= [rand*1000, rand*1000, 1000*1.1, 850*1.1];
P.font_size = 14;
P.axes_properties = {'FontSize', P.font_size, ...
                     'Color', 'none', ...
                     'TickDir', 'Out',...
                     'Nextplot', 'add', ...
                     'LineWidth', 1,...
                     'XColor', [0,0,0], ...
                     'YColor', [0,0,0], ...
                     'LabelFontSizeMultiplier', 1,...
                     'ActivePositionProperty','Position',...
                     'box','off'};  
P.panel_label_font_size = P.font_size * 1.5;
P.panel_label_pos = [0.1, 0.9, 0.1, 0.1];
P.panel_labels = char(65:90);
               