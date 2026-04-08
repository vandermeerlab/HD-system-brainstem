function [] = plotSaccadesAndAHVsinglePlot2(savedestination, varargin)

% 4/2021. JJS.
% Make a single plot that has saccade peths (both directions) and AHV tuning curve for a single neuron. Save all plots as images into a folder.
% This version includes a simple line for the PETHs (no errorbars) and the raw FR data for the tuning curve (more like a scatterplot).

% Inputs:
%   temporaldata: nCell x nBins double with firing rate data for temporal saccades  (from makeSaccadeHeatPlot.m)
%   nasaldata: nCell x nBins double with firing rate data for nasal saccades
%   *TC_data: not an explicit input, but this is retrieved from each session folder. Looks for file 'SSN-FRxAHV.mat'. Created by saveAHV_FR_scatterplot.m
%   binCenters: the centers of the bins used to calculate the saccade PETHs
%   cellname:  the name of each neuron for the saccade PETHs, in the order in which they were calculated (i.e. first folder in the directory)
%   cellID: the order of the saccade peth cells, relative to the order for all cells
%   nametouse: the name to use for saving each image file
%   savedestination: where the image files will be save to

occthresh = 0.5; % threshold number of seconds occupancy for including in tuning curve
formatSpec = '%.2f';
LineWidth = 5;
doSave = 1;
doPause = 1;
FontSize = 12;
process_varargin(varargin);
% formatSpec = '%.2f';

fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN;
    disp(SSN)
    
    if exist(strcat(SSN, '-VT1_proc.mat'))
        [FRxBinT, FRxBinN, FRxBinTsmooth, FRxBinNsmooth, FRxBinTnorm, FRxBinNnorm, TnormSmooth, NnormSmooth, outputIT, binCenters, cfg, cellID, cellname] = makeSaccadeHeatPlot;
        temporaldata = FRxBinTsmooth;
        nasaldata = FRxBinNsmooth;
        
        fc = FindFile('*FRxAHV.mat');    % this is pre-computed AHV tuning curve-type info. Should have variables AHV_tsd and AHV_F
        load(fc);
        [S] = LoadSpikesJeff;
        
        % get AHV Tuning Curve
        cfg_AHV = [];
        cfg_AHV.subsample_factor = 10;
        [AHV_tsd, tc_out] = AHV_tuning(S, cfg_AHV);
        AHV_dt = median(diff(AHV_tsd.tvec));
        
        for iPlot = 1:length(S.t)
            clf
            t = tiledlayout(1,1);
            ax1 = axes(t);
            plot(ax1, binCenters, nasaldata(iPlot,:), 'Color', [1 .6 .2], 'LineWidth', LineWidth);
            hold on
            plot(ax1, binCenters, temporaldata(iPlot,:), 'Color', 'r', 'LineWidth', LineWidth);
            legend('Nasal (CW)', 'Temporal (CCW)', 'Location', 'NorthWest')
            xlabel('Time peri Saccade (sec)')
            ylabel('Saccade Firing Rate (Hz)')
            ax1.XColor = 'r';
            ax1.YColor = 'r';
            ax2 = axes(t);
            %         plot(ax2, TC_data.bins(1,:),smoothdata(TC_data.tc(cellID(iPlot),:)), 'Color', 'k');
            plot(ax2, AHV_tsd.data, AHV_F(iPlot,:), '.', 'MarkerSize', .5, 'color', [.8 .8 .8]);
            hold on
            plot(ax2, tc_out.usr.binCenters(tc_out.occ_hist>occthresh), smoothdata(tc_out.tc(iPlot,(tc_out.occ_hist>occthresh))), 'k', 'LineWidth', 3);     
            c = axis;
            axis([-200 200 c(3) c(4)])
            ax2.XAxisLocation = 'top';
            ax2.YAxisLocation = 'right';
            ax2.Color = 'none';
            ax1.Box = 'off';
            ax2.Box = 'off';
            %     xlabel('AHV (deg./sec)')
            ylabel('AHV Firing Rate (Hz)')
            title(S.label{iPlot})
            
            h = get(gca, 'XLim');
            text(.75*h(1), 10, 'CW', 'FontSize', 12)
            text(.5*h(2), 10, 'CCW', 'FontSize', 12)
            
            %             text(.1,.85,num2str(zVal_T(iPlot),formatSpec), 'Units', 'Normalized', 'Color', 'r', 'FontSize', FontSize);
            %             text(.1,.65,num2str(zVal_N(iPlot),formatSpec), 'Units', 'Normalized', 'Color', [1 .6 .2], 'FontSize', FontSize);
            disp('paused. waiting for user')
            if doPause == 1
                pause
            end
            if doSave == 1
                baseFileName = S.label{iPlot};
                fname = fullfile(savedestination, strcat(baseFileName, '.eps'));
                saveas(gcf, fname)
                disp('saving')
            end
        end
    end
end
%
% M104-2021-01-24-FRxAHV.mat