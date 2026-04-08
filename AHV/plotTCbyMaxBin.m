[M, I] = max(TC_all.tc');

[B, isort] = sort(I);

pcolor(TC_all.bins(1,:), 1:length(TC_norm.tc), TC_norm.tc(isort, :))

cmap = jet(256);
% Make lowest one black
cmap(1,:) = 0;
colormap(cmap);