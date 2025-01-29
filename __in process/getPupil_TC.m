function [tc_pupil] = getPupil_TC(cfg_in, sd)
% JJS. 2023-04-28.
% Calculate and plot the tuning curve for eye position
% input:   sd - session data structure with spike trains S 
% output:  pupilTC - a 1 x nCell array of structures, each with the fields tc, occ_hist, spk_hist, usr, cfg 

cfg_def.doPlot = 1;
cfg_def.insetText = 18;
cfg_def.FontSize = 16;
cfg_def.occthresh = 2;   % occupancy threshold 
cfg_out = ProcessConfig2(cfg_def, cfg_in);
SSN = HD_GetSSN; disp(SSN);

try
    load(FindFile('*saccades-edited.mat'), 'tsdH') % tsdH is the pupil position variable
catch
    error('cannot find saccade data')
end

tsdH.data = tsdH.data'; % change the shape so that it is a "well-formed tsd" for tuning curves

for iCell = 1:length(sd.S.t)
    Stouse.t{1} = sd.S.t{iCell};
    cfg_tc = [];
%     cfg_tc.nBins = 50;
    cfg_tc.binEdges = {linspace(-60, 60, 101)};
    cfg_tc.occ_dt = median(diff(tsdH.tvec));
    cfg_tc.minOcc = 10;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
    tc_pupil{iCell} = TuningCurves(cfg_tc, Stouse, tsdH);
    
    if cfg_out.doPlot
        % calculate Q matrix
        tsdH_dt = median(diff(tsdH.tvec));
        cfg_Q = [];
        cfg_Q.smooth = 'gauss';
        cfg_Q.gausswin_sd = 0.05;
        cfg_Q.dt = tsdH_dt;
        cfg_Q.tvec_edges = tsdH.tvec(1):tsdH_dt:tsdH.tvec(end);
        F = MakeQfromS(cfg_Q, Stouse); % convert to FR
        F.data = F.data ./ cfg_Q.dt;
        
        % find FR corresponding to each AHV sample
        F_idx = nearest_idx3(tsdH.tvec, F.tvec);
        tsdH_F = F.data(:,F_idx);
        plot(tsdH.data, tsdH_F(1,:), '.', 'MarkerSize', .5, 'color', [.8 .8 .8]); hold on
        ylabel('Firing Rate (Hz)')
        
        plot(tc_pupil{iCell}.binCenters(tc_pupil{iCell}.occ_hist > cfg_out.occthresh), smoothdata(tc_pupil{iCell}.tc(1,(tc_pupil{iCell}.occ_hist > cfg_out.occthresh))), 'k', 'LineWidth', 3);
        set(gca, 'FontSize', cfg_def.FontSize)
        title('Pupil Position (pixels)')
        
        % axis tight
        c = axis;
        axis([-45 45 c(3) c(4)])
        line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
        text(-30, c(4)/2, 'nasal', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.15 .85 0])
        text(5, c(4)/2, 'temporal', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.55 .85 0])
        p.XAxisLocation = 'top';
    end
    disp('press any key to continue')
    pause
    clf
end
