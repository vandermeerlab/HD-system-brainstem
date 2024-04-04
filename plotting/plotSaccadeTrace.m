function plotSaccadeTrace(cfg_in)
% 2021-11. JJS.
% This function calculates saccades times from the eye position trace, similar to processPupilData2.feedback. Here, the threshold to use is manually adjusted,
%      and the trace is scrolled through to check/add/remove indivudal saccades that automatic method may have missed.
% This version accepts all potential saccades in the first pass, and sorts them (temporal vs. nasal) afterward, according to sign (+ = temp., - = nasal).
% This version is able to append new saccade entries/removals to the existing -saccade.m file data.

% INPUTS:
%           cfg_in: define variables such as the number of pixels
% OUTPUTS:
%           m.temporalSaccades:   timestamps for temporal saccades
%           m.nasalSaccades:      timestamps for nasal saccades
%           m.combinedSaccades:   both timestamps combined, sorted
%           m.tsdH:               tsd of horizontal pupil position
%           m.tsdV:               tsd of vertical pupil position
%           m.diffH:              tsd of horizontal pupil velocity  (*figure out units here)
%           m.diffV:              tsd of vertical pupil velocity
%           m.XS, YS:             x and y values for manually added TEMPORAL saccades
%           m.XN, YN:             x and y values for manually added NASAL saccades
%           m.temporalAmplitdues: amplitude of temporal saccades
%           m.nasalAmplitudes:    amplitude of nasal saccades
%           m.cfg:                record of parameters used for analysis, like thresholds

% 2024-03-01. JJS. Changed plot order and added figure legend.
SSN = HD_GetSSN; disp(SSN);

if exist(strcat(SSN, '-VT1.smi'), 'file') == 2
    SSN = HD_GetSSN;
    cfg_def = [];
    cfg_def.FontSize = 20;
    cfg_def.threshAdj  = 4;  % how many timesteps around a saccade that are disqualified from being considered as subsequent saccades. Saccades usually have a rebound that can hit the other threshold.
    cfg_def.threshT = 12;  % positive displacement in image pixel space. TEMPORAL saccades.
    cfg_def.threshN = -12; % negative displacement in image pixel space. NASAL saccades.
    
    cfg_def.scalingfactor = 1;  % for shrinking the pupil trace so its the same height as diffH
    cfg_def.artifactThresh = 4;  % units of pixels
    cfg_def.doPlotThresholds = 1;
    cfg_def.doPlotEverything = 1;
    cfg_def.LineWidth = 3;
    cfg_out = ProcessConfig2(cfg_def,cfg_in);  % not sure if there is a function difference btwn ver2 and ver1
    
    if exist(strcat(SSN, '-saccades-edited.mat'), 'file') == 2
        load(strcat(SSN, '-saccades-edited.mat'))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Plot the data and manually inspect
        
        hold on
        plot(diffH.tvec, diffH.data)
        plot(diffV.tvec, diffV.data, 'm')
        plot(tsdH.tvec, tsdH.data, 'Color', [.301 .745 .933], 'LineStyle', '--')
        %     plot(tsdH.tvec, tsdH.data, 'Color', 'k', 'LineStyle', '--')
        plot(amp1tsd.tvec, amp1tsd.data, 'k', 'LineWidth', cfg_out.LineWidth)
        plot(amp2tsd.tvec, amp2tsd.data, 'Color', [.85 .325 .098], 'LineWidth', cfg_out.LineWidth)
        xlabel('Time (sec)', 'FontSize', cfg_out.FontSize)
        ylabel('diff pupil pos', 'FontSize', cfg_out.FontSize)
%         title(sd.fc)
        line([tstart tend], [cfg_out.threshT cfg_out.threshT], 'Color', 'r')
        line([tstart tend], [cfg_out.threshN cfg_out.threshN], 'Color', 'g')
        line([tstart tend], [cfg_out.artifactThresh cfg_out.artifactThresh], 'Color', 'k', 'LineStyle', '--')
        line([tstart tend], [-cfg_out.artifactThresh -cfg_out.artifactThresh], 'Color', 'k', 'LineStyle', '--')
        plot(temporalSaccades, temporalAmplitudes, 'r.', 'MarkerSize', 25)
        plot(nasalSaccades, nasalAmplitudes, 'g.', 'MarkerSize', 25)
        set(gca, 'FontSize', cfg_out.FontSize)
        % legend('horiz eye vel.', 'vertical eye vel.', 'horizontal eye position', 'filtered vert. vel. 10-15 Hz', 'filtered horiz. vel. 10-15 Hz', '', '', '')
        yyaxis right
        plot(AHV_tsd.tvec, AHV_tsd.data, 'Color', [.75 .75 0])
        ylabel('horizontal pupil position', 'FontSize', cfg_out.FontSize)
        yyaxis left
    else
        disp('no saccade file detected')
    end
    
end