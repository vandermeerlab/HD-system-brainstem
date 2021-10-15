function [] = plotSaccadesAndAHVsinglePlot(temporaldata, nasaldata, TC_data, binCenters, cellname, cellCounterIndex, varargin)

% 4/2021. JJS.
% Make a single plot that has saccade peths (both directions) and AHV tuning curve for a single neuron. Save all plots as images into a folder.
% This version includes a simple line for the PETHs (no errorbars) and a simple line for the tuning curve (no errorbars). Does not include hemisphere info or zscores.

% savedestination = 'D:\Jeff\U01\analysis\figs\saccade PETH examples\saccade peth w AHV folder';
doSave = 1;
FontSize = 12;
process_varargin(varargin);
numCells = size(temporaldata,1);
formatSpec = '%.2f';

for iPlot = 1:numCells
    %     SSN = HD_GetSSN;
    %     disp(SSN)
    clf
    t = tiledlayout(1,1);
    ax1 = axes(t);
    plot(ax1, binCenters, temporaldata(iPlot,:), 'Color', 'r');
    hold
    plot(ax1, binCenters, nasaldata(iPlot,:), 'Color', [1 .6 .2]);
    legend('Temporal', 'Nasal', 'Location', 'NorthWest')
    xlabel('Time peri Saccade (sec)')
    ylabel('Saccade Firing Rate (Hz)')
    ax1.XColor = 'r';
    ax1.YColor = 'r';
    ax2 = axes(t);
    plot(ax2, TC_data.bins(1,:),smoothdata(TC_data.tc(cellCounterIndex(iPlot),:)), 'Color', 'k');
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right';
    ax2.Color = 'none';
    ax1.Box = 'off';
    ax2.Box = 'off';
%     xlabel('AHV (deg./sec)')
    ylabel('AHV Firing Rate (Hz)')
    title(cellname{iPlot})
    disp('paused. waiting for user')
    pause
    if doSave == 1
        saveas(gcf, num2str(iPlot), 'jpg')
        disp('saving')
    end
end