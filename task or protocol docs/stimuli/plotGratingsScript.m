objLeft = GratingStim([150 150], 'k', 50, 180, [.1 .3], 1, 0);
objRight = GratingStim([150 150], 'k', 50, 0, [.1 .3], 1, 0);

fullfig;
% pause(0.00001);
% frame_h = get(handle(gcf),'JavaFrame');
% set(frame_h,'Maximized',1);
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
plot(objLeft);


fullfig;
plot(objRight);

fullfig;
plot(objLeft);

fullfig;
plot(objRight);