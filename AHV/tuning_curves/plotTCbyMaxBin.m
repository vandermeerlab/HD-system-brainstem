[M, I] = max(TC_norm.tc(:,:),2);

[B, isort] = sort(I);

pcolor(TC_all.bins(1,:), 1:length(TC_norm.tc), TC_norm.tc(isort, :))