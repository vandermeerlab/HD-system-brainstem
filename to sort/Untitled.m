clf
for iNeuron = 1:104;
    plot(TC_out{iNeuron}.tc)
    title(num2str(iNeuron))
    disp('press any key')
    pause
end


box off
set(gca, 'FontSize', 18)
set(gca, 'TickDir', 'out');
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Color', 'k');


box off
set(gca, 'FontSize', 18)
set(gca, 'TickDir', 'out');
axis([-210 210 0 175]);
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-', 'Color', 'k');
ylabel('Firing Rate (Hz)')
xlabel('AHV (deg/s)')


