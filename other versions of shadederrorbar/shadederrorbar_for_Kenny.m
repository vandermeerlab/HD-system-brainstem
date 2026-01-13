
% brainstem code repo
clf
shadedErrorBar2(7:17, nanmean(AveMeanSpeedOTtouse_last(:,end-10:end)), nanstderr(AveMeanSpeedOTtouse_last(:,end-10:end))', {'r','markerfacecolor','r'});
hold on
shadedErrorBar2(14:17, nanmean(AveMeanSpeedCTtouse_last(:,end-3:end)), nanstderr(AveMeanSpeedCTtouse_last(:,end-3:end))', {'b','markerfacecolor','b'});


% plus maze code repo
% clf
% shadedErrorBar(14:17, nanmean(AveMeanSpeedCTtouse_last(:,end-3:end)), nanstderr(AveMeanSpeedCTtouse_last(:,end-3:end))', 'lineProps', 'b'); hold on
% shadedErrorBar(7:17, nanmean(AveMeanSpeedOTtouse_last(:,end-10:end)), nanstderr(AveMeanSpeedOTtouse_last(:,end-10:end))', 'lineProps', 'r')

clf
shadedErrorBar(8:11, nanmean(AveMeanSpeedCTtouse_last(:,end-3:end)), nanstderr(AveMeanSpeedCTtouse_last(:,end-3:end))', 'lineProps', 'b'); hold on
shadedErrorBar(1:11, nanmean(AveMeanSpeedOTtouse_last(:,end-10:end)), nanstderr(AveMeanSpeedOTtouse_last(:,end-10:end))', 'lineProps', 'r')

legend('CT', 'OT')
set(gca, 'XTick', [1:17])
set(gca, 'XTickLabel', '')
xlabel('Training Day')
ylabel('Speed (cm/s)')
set(gca, 'FontSize', 16)
axis([6 18 0 38]);
xlabel('Days Before Training Day')
set(gca, 'XTick', [11 10 9 8 7 6 5 4 3 2 1])