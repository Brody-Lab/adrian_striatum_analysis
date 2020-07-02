function [rat,session_date,data_fldr,recording_name] = extract_info_from_npx_path(path)
    session_date = regexprep(path,'.*(20[0-9][0-9][_-][0-9][0-9][_-][0-9][0-9]).*','$1');
    if strcmp(filesep,'\')
        rat = regexprep(path,'.*\\([A-Z][0-9][0-9][0-9])[\\_].*','$1');
    else
        rat = regexprep(path,'.*/([A-Z][0-9][0-9][0-9])[/_].*','$1');        
    end
    data_fldr = path(1:(strfind(path,rat)-1));
    [~,name] = fileparts(path);
    recording_name = strrep(name, '_g0', '');
end