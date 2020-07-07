function filename = fixfilesep(filename)

if strncmpi(filesep,'/',1)
    hits=strfind(filename,'\');
    filename(hits)=filesep;
else
    hits=strfind(filename,'/');
    filename(hits)=filesep;
end