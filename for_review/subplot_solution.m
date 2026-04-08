x1 = 0:0.1:40;
y1 = 4.*cos(x1)./(x1+2);
x2 = 1:0.2:20;
y2 = x2.^2./x2.^3;f = figure;
ax1 = subplot(2,1,2);
ax2 = copyobj(ax1, 1);
plot(ax1,x1,y1,'-r')
ax1.XColor = 'r';
ax1.YColor = 'r';plot(ax2,x2,y2,'-k')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';