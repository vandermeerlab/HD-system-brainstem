function [data_out, data_outR, nSpikesRemoved] = getSlowPhaseData(cfg_in, sd, iCell)
% function [data_out, nSpikesRemoved] = getSlowPhaseData(cfg_in, sd, iCell)
%
%   Inputs
%           iCell - which cell number to plot. leave empty [] if not calculating neural data
%
% use saccade .mat file to restrict data to slow phase only
% useful for later tuning curve and GLM analysis
%
% 2024-04-18. MvdM, based on JJS getSlowPhaseTC3
% 2024-05-29. JJS. Added conditions to skip the resampling of neural data.

cfg_def = [];
cfg_def.FontSize = 20;
cfg_def.saccade_pre = 0.2;  % how many seconds to cut out before the saccade.
cfg_def.saccade_post = 0.1; % how many seconds to cut out after the saccade.
cfg_def.doPlot = 1;
cfg_def.markerSize = 3;
cfg_def.plotWhichPhase = 'slow'; % 'slow', 'both'
cfg_def.medfilt_window_size = 11; % number of samples to use for slow phase velocity filter
cfg_def.smooth = 1;

cfg_Q = [];
cfg_Q.dt = 0.02; % binsize in s
cfg_Q.smooth = 'gauss'; % [], 'gauss'
cfg_Q.gausswin_size = 1; % gaussian window size in seconds; only used if cfg.smooth = 'gauss'
cfg_Q.gausswin_sd = 0.1; % SD for gaussian convolution in seconds
cfg_def.cfg_Q = cfg_Q;
nSpikesRemoved = 0;

cfg_master = ProcessConfig2(cfg_def, cfg_in);

if exist(strcat(sd.SSN, '-saccades-edited.mat'),'file') == 2
    load(strcat(sd.SSN, '-saccades-edited.mat')) % this will load the edited (curated) saccade timestamps and the pupil trace time series (position and velocity)
else
    error('saccade mat file not found')
end

if ~isempty(iCell)
    this_cell = SelectTS([], sd.S, iCell); % choose a single neuron to work with
    meanFRoverall = length(this_cell.t{1}) / sd.SessLength;
end

%% resample data to be on common timebase
new_tvec_edges = sd.AHV.tvec(1):cfg_master.cfg_Q.dt:sd.AHV.tvec(end);
new_tvec = new_tvec_edges(1:end-1) + (cfg_master.cfg_Q.dt / 2);

% compute firing rate
if ~isempty(iCell)
    cfg_master.cfg_Q.tvec_edges = new_tvec_edges;
    fr = MakeQfromS(cfg_master.cfg_Q, this_cell); fr.data = fr.data ./ cfg_master.cfg_Q.dt;
end

% ahv
cfg = []; cfg.method = 'rebinning';
ahv = resampleTSD(cfg, sd.AHV, new_tvec);

% eye position
horiz_eye_pos = resampleTSD(cfg, sd.tsdH, new_tvec);

% eye velocity
%horiz_eye_vel = dxdt(sd.tsdH.tvec, sd.tsdH.data, 'window', 0.2, 'postSmoothing', 0.1); % this would be nice, why doesn't it work?
sd.horiz_eye_vel = cat(2, 0, diff(sd.tsdH.data)); sd.horiz_eye_vel = sd.horiz_eye_vel ./ median(diff(sd.tsdH.tvec));
sd.horiz_eye_vel_smooth = medfilt1(sd.horiz_eye_vel, cfg_master.medfilt_window_size);
sd.horiz_eye_vel = tsd(sd.tsdH.tvec, sd.horiz_eye_vel);
sd.horiz_eye_vel_smooth = tsd(sd.tsdH.tvec, sd.horiz_eye_vel_smooth);

if cfg_master.doPlot
    clf;
    subplot(221)
    displayFactor = 10;
    plot(sd.tsdH.tvec, sd.tsdH.data .* displayFactor); hold on  % pupil position
    plot(sd.horiz_eye_vel,'r');  % pupil velocity
    plot(sd.horiz_eye_vel_smooth, 'g', 'LineWidth', 2) % smoothed pupil velocity
    set(gca, 'FontSize', cfg_master.FontSize)
    plot(sd.horiz_eye_vel,'r.')
    plot(sd.tsdH.tvec, sd.tsdH.data .* displayFactor, '.b')
    plot(sd.horiz_eye_vel_smooth, 'g', 'LineWidth', 2) % smoothed pupil velocity. Replotted here so the green line is in front.
    legend('pupil position', 'velocity', 'smoothed velocity')
end

cfg = []; cfg.method = 'rebinning';
horiz_eye_vel = resampleTSD(cfg, sd.horiz_eye_vel, new_tvec);
horiz_eye_vel_smooth = resampleTSD(cfg, sd.horiz_eye_vel_smooth, new_tvec);

% head direction
cfg = []; cfg.method = 'rebinning';
head_direction = resampleTSD(cfg, sd.orientation, new_tvec);

%% Limit the data to everything but the quick phase periods
timestouse = ~isnan(combinedSaccades); % create a logical with not NaN saccade times
combinedSaccadesToUse = combinedSaccades(timestouse); % select only non-NaN saccade values to use

pre = combinedSaccadesToUse - cfg_master.saccade_pre;
post = combinedSaccadesToUse + cfg_master.saccade_post;

if ~isempty(iCell)
    [this_cellR, keep] = antirestrict(this_cell, pre, post);
    nSpikes = length(this_cell.t{1});
    fprintf('Total SPIKES: %d\n', nSpikes);
    nSpikesRemoved = sum((keep{1} == 0));
    fprintf('num SPIKES removed: %d (fraction %.2f)\n', nSpikesRemoved, nSpikesRemoved/nSpikes);
    frR = antirestrict(fr, pre, post);
end

%% restrict everything else and set up data_out
ahvR = antirestrict(ahv, pre, post);
horiz_eye_posR = antirestrict(horiz_eye_pos, pre, post);
horiz_eye_velR = antirestrict(horiz_eye_vel, pre, post);
horiz_eye_vel_smoothR = antirestrict(horiz_eye_vel_smooth, pre, post);
head_directionR = antirestrict(head_direction, pre, post);

if cfg_master.doPlot
    subplot(222)
    plot(ahvR.data, horiz_eye_vel_smoothR.data,  '.b', 'MarkerSize', cfg_master.markerSize);
    xlabel('AHV'); ylabel('horiz eye vel. (smoothed)')
    set(gca, 'FontSize', cfg_master.FontSize)
    lsline
end

% restricted and resampled data
data_outR.horiz_eye_pos = horiz_eye_posR; data_out.horiz_eye_pos = horiz_eye_pos;
data_outR.horiz_eye_vel = horiz_eye_velR; data_out.horiz_eye_vel = horiz_eye_vel;
data_outR.horiz_eye_vel_smooth = horiz_eye_vel_smoothR; data_out.horiz_eye_vel_smooth = horiz_eye_vel_smooth;
data_outR.AHV = ahvR; data_out.AHV = ahv;
data_outR.HD = head_directionR; data_out.HD = head_direction;

if ~isempty(iCell)
    data_outR.fr = frR; data_out.fr = fr;
    data_outR.S = this_cellR; data_out.S = this_cell;
end
if cfg_master.doPlot
    
    cfg_master.LineWidth = 3;
    %     cfg_master.FontSize = 20;
    cfg_master.insetText = 18;
    cfg_master.tightX = .025;
    cfg_master.tightY = .025;
    cfg_master.smooth = 1;
    
    switch cfg_master.plotWhichPhase
        case 'slow'
            this_data = data_outR; % can use this to set to something else
            title_string = 'Slow Phase Only';
        case 'both'
            this_data = data_out;
            title_string = 'Slow + Quick Phase';
    end
    
    %% get AHV Tuning Curve
    cfg_tcAHV = [];
    cfg_tcAHV.nBins = 100;
    cfg_tcAHV.binEdges = {linspace(-200, 200, 101)};
    cfg_tcAHV.minOcc = 5;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
    cfg_tcAHV.occ_dt = median(diff(ahvR.tvec));
    tc_outAHV = TuningCurves(cfg_tcAHV, this_data.S, this_data.AHV);
    
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    %% #1 AHV tuning curve
    % calculate raw firing rates
    %     figure; clf;
    hold on;
    subplot(2,2,3)
    %     subtightplot(4,2,1, [cfg_master.tightX cfg_master.tightY]);
    
    % NOTE it is critical to not recompute firing rates here (as in
    % previous versions!)
    F_idx = nearest_idx3(this_data.AHV.tvec, this_data.fr.tvec);
    AHV_F = this_data.fr.data(:,F_idx);
    
    ymax = max(AHV_F);
    set(gca, 'TickDir', 'out', 'FontSize', cfg_master.FontSize)
    plot(this_data.AHV.data, AHV_F, '.b', 'MarkerSize', cfg_master.markerSize); hold on
    set(gca, 'Ylim', [0 ymax], 'FontSize', cfg_master.FontSize)
    
    % add the Tuning Curve
    if cfg_master.smooth
        plot(tc_outAHV.binCenters, smoothdata(tc_outAHV.tc), 'LineWidth', cfg_master.LineWidth, 'Color', 'k');
    else
        plot(tc_outAHV.binCenters, tc_outAHV.tc, 'LineWidth', 3, 'Color', 'k');
    end
    xlabel('AHV (deg/s)', 'FontSize', cfg_master.FontSize)
    ylabel('FR (Hz)', 'FontSize', cfg_master.FontSize)
    set(groot, 'DefaultLegendInterpreter', 'none')
    title(title_string)
    text(NaN, NaN, 'CW', 'FontSize', cfg_master.insetText, 'Units', 'normalized', 'Position', [.25 .85 0])
    text(NaN, NaN, 'CCW', 'FontSize', cfg_master.insetText, 'Units', 'normalized', 'Position', [.6 .85 0])
    c = axis;
    line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
    
    %% #2 Eye Velocity Tuning Curve
    cfg_V = [];
    cfg_V.nBins = 100;
    cfg_V.binEdges = {linspace(-200, 200, 101)};
    cfg_V.occ_dt = median(diff(this_data.horiz_eye_vel_smooth.tvec));
    cfg_V.minOcc = 5;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
    tc_velEV = TuningCurves(cfg_V, this_data.S, this_data.horiz_eye_vel_smooth);
    
    % Calculate pupil velocity firing rates
    subplot(2,2,4)
    %     subtightplot(4,2,2, [cfg_master.tightX cfg_master.tightY]);
    hold on
    % find FR corresponding to each pupil position sample
    F_idxEV = nearest_idx3(this_data.horiz_eye_vel_smooth.tvec, this_data.fr.tvec);
    EV_F = this_data.fr.data(:,F_idxEV);
    ymaxEV = max(EV_F);
    plot(this_data.horiz_eye_vel_smooth.data, EV_F, '.b', 'MarkerSize', cfg_master.markerSize);   %  'color', [.8 .8 .8]
    
    set(gca, 'FontSize', cfg_master.FontSize)
    % title('Pupil Position (pixels)')
    % axis tight
    c = axis;
    axis([-45 45 c(3) c(4)])
    line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
    text(-30, c(4)/2, 'nasal', 'FontSize', cfg_master.insetText, 'Units', 'normalized', 'Position', [.15 .85 0])
    text(30, c(4)/2, 'temporal', 'FontSize', cfg_master.insetText, 'Units', 'normalized', 'Position', [.75 .85 0])
    xlabel('Eye Velocity (pixels/s)', 'FontSize', cfg_master.FontSize)
    ylabel('FR (Hz)', 'FontSize', cfg_master.FontSize)
    set(gca, 'FontSize', cfg_master.FontSize)
    title(sd.fn{iCell,1})
    
    if cfg_master.smooth
        plot(tc_velEV.binCenters, smoothdata(tc_velEV.tc), 'LineWidth', 3, 'Color', 'k');
    else
        plot(tc_velEV.binCenters, tc_velEV.tc, 'LineWidth', 3, 'Color', 'k');
    end
    c = axis;
    axis([-150 150 c(3) max(this_data.fr.data)]); c = axis;
    line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
    
end