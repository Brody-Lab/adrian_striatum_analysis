function region_flag = is_in_region(Cells,cellno,regions) 
    if isempty(cellno)
        cellno=1:sum(Cells.n_clusters);
    end
    ncells = numel(cellno);
    region_flag = false(size(cellno));
    if numel(Cells.n_clusters)>1
        pen_idx = get_pen_idx(cellno,Cells.n_clusters);
    else
        pen_idx=ones(ncells,1);
    end
    for i=1:ncells
        if Cells.regions(cellno(i))
            names = Cells.penetration(pen_idx(i)).regions(Cells.regions(cellno(i))).name;
            if any(ismember(regions,names))
                region_flag(i)=true;
            end
        end
    end
end

function pen_idx = get_pen_idx(cellno,cell_count)
    ncells=numel(cellno);
    pen_idx = ones(size(cellno));
    if numel(cell_count)>2
        error('Case of 3 joined cells file not implemented yet.');
    end
    for i=1:ncells
        if cellno(i)>cell_count(1)
            pen_idx(i) = 2;
        end
    end
end