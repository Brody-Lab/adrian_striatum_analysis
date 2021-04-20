function plot_cell_position()
    P=get_parameters;
    T=load_cells_table();
    T.ML(T.hemisphere=="left") = -T.ML(T.hemisphere=="left");
    scatter3(T.ML,T.AP,T.DV,'.');
    ylabel('AP (mm anterior to bregma)');
    zlabel('DV (mm below brain surface');
    xlabel('ML (mm right of midline');
    set(gca,'zdir','reverse',P.axes_properties{:},'CameraPosition',[39 -36 -26]);grid on;
end