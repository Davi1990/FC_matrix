function view_FC(zFC,parcels_list)

[~, cont, ~, COMM_list_red] = make_Gordon_parcels_table(parcels_list);
tick_label = [0 cumsum(cont(1:end-1))] + round(cont./2);

figure, imagesc(zFC,[-0.4, 1.2]),
title('fMRI (zFisher-transformed) FC map')
colormap jet
colorbar
axis square
set(gca,'XTick',tick_label, 'XTickLabel', COMM_list_red,'XTickLabelRotation',45)
set(gca,'YTick',tick_label, 'YTickLabel', COMM_list_red)
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'Color',[1 1 1]);
