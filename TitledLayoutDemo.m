


t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1, binCenters, temporaldata(iPlot,:), 'Color', 'r');
hold
plot(ax1, binCenters, nasaldata(iPlot,:), 'Color', 'r');
ax1.XColor = 'r';
ax1.YColor = 'r';
ax2 = axes(t);
plot(ax2, TC_all.bins(1,:),smoothdata(TC_all.tc(iPlot,:)), 'Color', 'k');
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
