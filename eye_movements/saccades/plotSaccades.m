function [] = plotSaccades(AHV_tsd, amp1tsd, amp2tsd, tsdH, diffH, diffV, temporalSaccades, temporalAmplitudes, nasalSaccades, nasalAmplitudes)

LineWidth = 3;
FontSize = 12;
tstart = diffH.tvec(1);
tend = diffH.tvec(end);
threshT = 10;
threshN = 10; 
artifactThresh = 4; 
%% Plot the data and manually inspect
SSN = HD_GetSSN; disp(SSN);
clf;
hold on
plot(diffH.tvec, diffH.data)
plot(diffV.tvec, diffV.data, 'm')
plot(tsdH.tvec, tsdH.data, 'Color', [.301 .745 .933], 'LineStyle', '--')
plot(amp1tsd.tvec, amp1tsd.data, 'k', 'LineWidth', LineWidth)
plot(amp2tsd.tvec, amp2tsd.data, 'Color', [.85 .325 .098], 'LineWidth', LineWidth)
xlabel('Time (sec)', 'FontSize', FontSize)
ylabel('diff pupil pos', 'FontSize', FontSize)
title(SSN)
line([tstart tend], [threshT threshT], 'Color', 'r')
line([tstart tend], [threshN threshN], 'Color', 'g')
line([tstart tend], [artifactThresh artifactThresh], 'Color', 'k', 'LineStyle', '--')
line([tstart tend], [-artifactThresh -artifactThresh], 'Color', 'k', 'LineStyle', '--')
plot(temporalSaccades, temporalAmplitudes, 'r.', 'MarkerSize', 25)
plot(nasalSaccades, nasalAmplitudes, 'g.', 'MarkerSize', 25)
set(gca, 'FontSize', FontSize)
% legend('horiz eye vel.', 'vertical eye vel.', 'horizontal eye position', 'filtered vert. vel. 10-15 Hz', 'filtered horiz. vel. 10-15 Hz', '', '', '')
yyaxis right
plot(AHV_tsd.tvec, AHV_tsd.data, 'Color', [.75 .75 0])
% c = axis;
% line([c(1) c(2)], [0 0], 'Color', [.75 .75 0], 'LineStyle', '--', 'LineWidth', 3)        % plotting another line makes it glitchy with the
ylabel('horizontal pupil position', 'FontSize', FontSize)
yyaxis left