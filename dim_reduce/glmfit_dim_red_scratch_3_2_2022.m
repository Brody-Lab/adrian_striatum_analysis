ap_group = equalize_ap_groups(fits_table.ap_group);


[coeff, score,latent, tsquared, explained] = pca(bs,'NumComponents',3);

figure;
colors=parula(4);
for i=1:4
    idx= ap_group==i;
    scatter((coeff(idx,1)),(coeff(idx,2)),'CData',colors(i,:));hold on;
end