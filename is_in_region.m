function region_flag = is_in_region(Cells,cellno,regions) 
    if isempty(cellno)
        cellno=1:sum(Cells.n_clusters);
    end
    region_flag = ismember(Cells.region(cellno),regions);    
end