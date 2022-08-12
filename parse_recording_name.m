function [rat,date,suffix] = parse_recording_name(recording_name)
    pattern='([A-Z]\d{3})_(\d{4}_\d{2}_\d{2})(.*)';
    rat = regexprep(recording_name,pattern,'$1');
    date = regexprep(recording_name,pattern,'$2');
    suffix = regexprep(recording_name,pattern,'$3');
    if any(strcmp(rat,recording_name) | strcmp(date,recording_name) | strcmp(suffix,recording_name))
        error('%s does not match the expected pattern of a recording name: rat_yyyy_mm_dd...',recording_name);
    end
end
