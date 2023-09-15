P=get_parameters;

fid = fopen([P.repository_path,filesep,'ap_group_colors.json'],'w');
n=numel(P.ap_groups);

txt = jsonencode(struct('ap_group',mat2cell((1:n)',ones(1,n),1),'rgb_color',mat2cell(P.ap_group_colors,ones(1,n))),'PrettyPrint',true);
fprintf(fid,txt);

fclose(fid);