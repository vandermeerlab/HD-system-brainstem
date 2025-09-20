function fingerprintPlot(cfg_in, sd, iCell)
% 2023-07-19. JJS.
%   This function creates a subplot figure for each individual neuron within a session folder.
%   It calculates and displays all of the relevant features for each neuron. If data is absent or not applicable, those parts of the figure will be blank.

% Inputs:
%       cfg     - variable inputs. If cfg is blank, this function will use default parameters.
%       sd      - structure with session data for the given folder. Use sd = LoadSessionData(fd).
%       iCell   - which neuron to plot (goes in numerical order from TT1)
% Outputs:
%       optional write file

clf;    % clear previous figure, if present
tic;    % keep track of how long this function takes
%% Initialize Parameters
cfg_def = [];
cfg_def.tightX = .025;      % related to subtightplot.  X and Y = 'gap'- a two elements vector [vertical,horizontal] defining the gap between neighbouring axes.
cfg_def.tightY = .02;       % see comment above
cfg_def.occthresh = 0.5;    % Occupancy threhsold for including data in pupil position tuning curve
cfg_def.smallfont = 8;      % Font size for the saccade peth legends.
cfg_def.insetText = 18;     % Font size of the inset text ('CW', 'CCW', or 'nasal', 'temporal') for the first two subplots.
cfg_def.speedthresh = 0.3;  % Wheel velocity values that are in between +- cfg.speedthresh (cm/s) are not displayed in the wheel speed TC subplot. These are very low speed values where the mouse is essentially stationary. Including them introduces noise into the regression line.
cfg_def.ahv_thresh = 4;     % AHV values that are in between +- cfg.ahv_thresh will be excluded from WHEEL speed vs. AHV scatterplot. These are very low AHV values and introduce noise into the regression line.
cfg_def.doLaser = 1;        % Choose 0 to exclude laser events. Choose 1 to inlude laser events. This is a temporary measure b/c there are bugs in some sessions for laser event time extraction
cfg_def.FontSize = 13;      % General font size
cfg_def.histXmin = 0.01;    % Axis size for the HISTISI plot
cfg_def.histXmax = 0.2;
cfg_def.LineWidth = 3;      % General line width for plotting

cfg = ProcessConfig2(cfg_def, cfg_in);

%% Is there Eyetracking data for this session?
if exist(strcat(sd.SSN, '-VT1_proc.mat'))
    eye = 1;
else
    eye = 0;
end
%% Isolate the spike train for the single neuron indicated in the function input. This is so that restrict() can operate within spikePETH()
spikefiles = FindFiles('*.t');
cfg.fc = {spikefiles{iCell}}; % limit to one neuron
% S = LoadSpikesJeff(cfg);  % using this loader to subtract out starttime from unix time   %% JJS. 2025-09-20. This is superceded by sd = LoadSession

%% Change the .t filename from having an underscore '_' to a dash '-', for ease of use later. Most filenames use a dash '-'.
[fc] = FindFiles(strcat(sd.SSN, '*.t'));
[~, b, ~] = fileparts(fc);
if iscell(b)
    cellID = b{iCell};
else
    cellID = b;
end
newID = cellID;
k = strfind(cellID, '_');
newID(k) = '-';

%% #1 Firing Rate x AHV (with Tuning Curve)
p = subtightplot(6,6,1, [cfg.tightX cfg.tightY]);
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = sd.AHV.dt;
cfg_Q.tvec_edges = sd.AHV.tvec(1): sd.AHV.dt: sd.AHV.tvec(end);
F = MakeQfromS(cfg_Q, sd.S); % convert to FR
F.data = F.data(iCell,:) ./ cfg_Q.dt;
F_idx = nearest_idx3(sd.AHV.tvec, F.tvec);
AHV_F = F.data(:,F_idx);
ymax = max(AHV_F);
set(gca, 'TickDir', 'out', 'FontSize', cfg.FontSize)
plot(sd.AHV.data, AHV_F, '.', 'MarkerSize', .5); hold on
set(gca, 'Ylim', [0 ymax], 'FontSize', cfg.FontSize)

% Add Tuning Curve
cfg_ahv.doPlot = 0;
[tc_out] = getAHV_TC(cfg_ahv, sd);
plot(tc_out.binCenters, smoothdata(tc_out.tc(iCell,:)), 'LineWidth', cfg.LineWidth, 'Color', 'k');
ylabel('FR (Hz)', 'FontSize', cfg.FontSize)
set(groot, 'DefaultLegendInterpreter', 'none')
title('AHV Tuning Curve')
text(NaN, NaN, 'CW', 'FontSize', cfg.insetText, 'Units', 'normalized', 'Position', [.25 .85 0])
text(NaN, NaN, 'CCW', 'FontSize', cfg.insetText, 'Units', 'normalized', 'Position', [.6 .85 0])
p.XAxisLocation = 'top';
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
axis tight
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #2 Firing Rate x Pupil Position (with Tuning Curve)
if eye
    
    p = subtightplot(6,6,2, [cfg.tightX cfg.tightY]); hold on
    % calculate Q matrix
    cfg_Q = [];
    cfg_Q.smooth = 'gauss';
    cfg_Q.gausswin_sd = 0.05;
    cfg_Q.dt = sd.tsdH_dt;
    cfg_Q.tvec_edges = sd.tsdH.tvec(1): sd.tsdH_dt: sd.tsdH.tvec(end);
    F = MakeQfromS(cfg_Q, sd.S); % convert to FR
    F.data = F.data(iCell,:) ./ cfg_Q.dt;
    
    % find FR corresponding to each pupil position sample
    F_idx = nearest_idx3(sd.tsdH.tvec, F.tvec);
    tsdH_F = F.data(:,F_idx);
    plot(sd.tsdH.data, tsdH_F(1,:), '.', 'MarkerSize', .5, 'color', [.8 .8 .8]); hold on
    
    tsdH.data = sd.tsdH.data'; % change the shape so that it is a "well-formed tsd" for tuning curves
    cfg_tc = [];
    cfg_tc.nBins = 50;
    cfg_tc.binEdges = {linspace(-60, 60, 101)};
    cfg_tc.occ_dt = median(diff(sd.tsdH.tvec));
    cfg_tc.minOcc = 10;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
    tc_pupil = TuningCurves(cfg_tc, sd.S, sd.tsdH);
    plot(tc_pupil.binCenters(tc_pupil.occ_hist > cfg.occthresh), smoothdata(tc_pupil.tc(1,(tc_pupil.occ_hist > cfg.occthresh))), 'k', 'LineWidth', 3);
    set(gca, 'FontSize', cfg.FontSize)
    title('Pupil Position (pixels)')
    % axis tight
    c = axis;
    axis([-45 45 c(3) c(4)])
    line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
    text(-30, c(4)/2, 'nasal', 'FontSize', cfg.insetText, 'Units', 'normalized', 'Position', [.15 .85 0])
    text(5, c(4)/2, 'temporal', 'FontSize', cfg.insetText, 'Units', 'normalized', 'Position', [.55 .85 0])
    p.XAxisLocation = 'top';
    title('Pupil Position (pixels)')
    clear F
end
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #3 acf
p = subtightplot(6,6,3, [cfg.tightX cfg.tightY]);
cfg_acf = [];
cfg_acf.binsize = 0.001;
cfg_acf.max_t = 0.5;
cfg_acf.smooth = 1; % set to 1 to compute ccf on SDF, 0 for raw spikes
cfg_acf.gauss_w = 1; % width of Gaussian convolution window (in s)
cfg_acf.gauss_sd = 0.005; % SD of Gaussian convolution window (in s)
set(gca, 'TickDir', 'out', 'FontSize', cfg.FontSize)
[acf, tvec] = ccf(cfg_acf, sd.S.t{1}, sd.S.t{1});
midpoint = ceil(length(acf)./2);
acf(midpoint) = 0;
plot(tvec, acf, 'LineWidth', 1);
% xlabel('Time (sec)', 'FontSize', FontSize)
set(gca, 'xtick', [-.5 -.25 0 .25 .5], 'FontSize', cfg.FontSize); grid on;
title('Acorr')
p.XAxisLocation = 'top';
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #4 acf zoomed in
p = subtightplot(6,6,4, [cfg.tightX cfg.tightY]); hold on
cfg_acf = [];
cfg_acf.binsize = 0.001;
cfg_acf.max_t = 0.05;
cfg_acf.smooth = 1; % set to 1 to compute ccf on SDF, 0 for raw spikes
cfg_acf.gauss_w = 1; % width of Gaussian convolution window (in s)
cfg_acf.gauss_sd = 0.005; % SD of Gaussian convolution window (in s)
set(gca, 'TickDir', 'out', 'FontSize', cfg.FontSize)
[acf, tvec] = ccf(cfg_acf, sd.S.t{1}, sd.S.t{1});
midpoint = ceil(length(acf)./2);
acf(midpoint) = 0;
plot(tvec, acf, 'LineWidth', 1);
% xlabel('Time (sec)', 'FontSize', FontSize)
set(gca, 'xtick', [-.05 0 .05], 'FontSize', cfg.FontSize); grid on;
title('Acorr')
p.XAxisLocation = 'top';
%---------------------------------------------------------------------------------------------------------------------------------------------------
%% #5   Wheel speed Tuning Curve
p = subtightplot(6,6,5, [cfg.tightX cfg.tightY]); hold on
ylabel('FR (Hz)', 'FontSize', cfg.FontSize)
title('Wheel Speed Tuning Curve', 'FontSize', cfg.FontSize)

WheelencoderCSC = strcat(sd.SSN, '-CSC34.Ncs');  % check to see if this session has running wheel data
if exist(WheelencoderCSC)
    % calculate Q matrix
    cfg_Q = [];
    cfg_Q.smooth = 'gauss';
    cfg_Q.gausswin_sd = 0.05;
    cfg_Q.dt = sd.speed.dt;
    cfg_Q.tvec_edges = sd.speed.tvec(1): sd.speed.dt: sd.speed.tvec(end);
    F = MakeQfromS(cfg_Q, sd.S); % convert to FR
    F.data = F.data(iCell,:) ./ cfg_Q.dt;
    
    % find FR corresponding to each AHV sample
    F_idx = nearest_idx3(sd.speed.tvec, F.tvec);
    tsdH_F = F.data(:,F_idx);
    yyaxis right
    
    % Plot the raw speed data
    z = sd.speed.data > cfg.speedthresh | sd.speed.data < -cfg.speedthresh;  % remove data at very low speeds
    plot(sd.speed.data(z), tsdH_F(1,z), '.', 'MarkerSize', .5, 'color', [.8 .8 .8]); hold on
    h = lsline;
    set(h(1), 'Color', 'k')
    set(h(1), 'LineWidth', 2)
    
    % Plot a regression line, overlaid
    speeddata = sd.speed.data(z)';
    speeddata(:,2) = ones(length(speeddata), 1);
    [~,~,~,~,stats] = regress(tsdH_F(1,z)', speeddata);
    c = axis;
    line([-cfg.speedthresh -cfg.speedthresh], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
    line([cfg.speedthresh cfg.speedthresh], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
    text(NaN, NaN, strcat('Rsq =', sprintf('%0.2f', stats(1))), 'FontSize', 12, 'Units', 'normalized', 'Position', [.55 .85 0])
    
    % Plot the tuning curve, overlaid
    cfg_tc = [];
    cfg_tc.nBins = 50;
    cfg_tc.binEdges = {linspace(-5, 30, 101)};
    cfg_tc.occ_dt = median(diff(sd.speed.tvec));
    cfg_tc.minOcc = 1;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
    tc_speed = TuningCurves(cfg_tc, sd.S, sd.speed);
    plot(tc_speed.binCenters(tc_speed.occ_hist > cfg.occthresh), smoothdata(tc_speed.tc(1,(tc_speed.occ_hist > cfg.occthresh))), 'k', 'LineWidth', 3);
    axis tight
    yyaxis left
    set(gca, 'YTick', [])
    yyaxis right
    grid on
    set(gca, 'FontSize', cfg.FontSize)
    p.XAxisLocation = 'top';
    axis tight
end

%% #6  WHEEL speed vs. AHV scatterplot
plot6 = subtightplot(6,6,6, [cfg.tightX cfg.tightY]);
xlabel('AHV', 'FontSize', cfg.FontSize)
ylabel('Wheel Speed', 'FontSize', cfg.FontSize)
title('Wheel Speed Tuning Curve', 'FontSize', cfg.FontSize)

if exist(WheelencoderCSC)              % see if this session has wheel speed data
    X = sd.AHV.data > cfg.ahv_thresh | sd.AHV.data < -cfg.ahv_thresh;  % remove the wheel data where AHV is close to zero (i.e., platform is stationary)
    Z_idx = nearest_idx3(sd.AHV.tvec(X), sd.speed.tvec);
    SPEED = tsd(sd.speed.tvec(Z_idx,:), sd.speed.data(:,Z_idx));
    plot(sd.AHV.data(X), SPEED.data, '.', 'MarkerSize', 1);
    axis tight
    xlabel('AHV', 'FontSize', cfg.FontSize)
    ylabel('Wheel Speed', 'FontSize', cfg.FontSize)
    % Plot a regression line and the threshold lines
    h = lsline;
    set(h(1), 'Color', 'k')
    set(h(1), 'LineWidth', 2)
    c = axis;
    line([c(1) c(2)], [0 0], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
    line([-cfg.ahv_thresh -cfg.ahv_thresh], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
    line([cfg.ahv_thresh cfg.ahv_thresh], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
    % She the r value as text
    SPEEDregress = SPEED.data';
    SPEEDregress(:,2) = ones(length(SPEED.data), 1);
    [~,~,~,~,stats] = regress(sd.AHV.data(X)', SPEEDregress);
    text(NaN, NaN, strcat('Rsq =', sprintf('%0.2f', stats(1))), 'FontSize', 12, 'Units', 'normalized', 'Position', [.05 .85 0])
    plot6.YAxisLocation = 'right';
    plot6.XAxisLocation = 'top';
    set(gca, 'FontSize', cfg.FontSize)
end
%% #7 MOVING SACCADE peth: wide
if eye
    subtightplot(6,6,7, [cfg.tightX cfg.tightY]); hold on
    cfg_1.doPlot = 0;
    cfg_1.window = [-2 2];
    cfg_1.dt = 0.05;
    [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, sd.S, sd.nasal_timestamps_MOVING);
    [mn1, edges] = histcounts(outputS_n, outputIT_n);
    plot(edges(1:end-1), mn1/cfg_1.dt/length(sd.nasal_timestamps_MOVING), 'LineWidth', cfg.LineWidth); % this is a hack. should replace binedgges with bincenters
    hold on
    [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, sd.S, sd.temporal_timestamps_MOVING);
    [mt1, edges] = histcounts(outputS_t, outputIT_t);
    plot(edges(1:end-1), mt1/cfg_1.dt/length(sd.temporal_timestamps_MOVING), 'LineWidth', cfg.LineWidth);
    c = axis;
    line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
    set(gca, 'FontSize', cfg.FontSize)
    legend('nasal', 'temporal', '', 'FontSize', cfg.smallfont)
    set(gca, 'XTick', [-2 2])
    ylabel('FR (Hz)', 'FontSize', cfg.FontSize)
    
    % text(NaN, NaN, sprintf('Lratio %.2f \n', CluSep.Lratio), 'FontSize', 10, 'Units', 'normalized', 'Position', [.2 .85 0])
    text(NaN, NaN, strcat('#temp.=', sprintf('%0.0f', sd.numTemporal_moving)), 'FontSize', 16, 'Units', 'normalized', 'Position', [.02 .95 0])
    text(NaN, NaN, strcat('#nasal=', sprintf('%0.0f', sd.numNasal_moving)), 'FontSize', 16, 'Units', 'normalized', 'Position', [.02 .80 0])
    
    %% #8 Stationary SACCADE peth: wide
    subtightplot(6,6,8, [cfg.tightX cfg.tightY]); hold on
    if ~isempty(sd.nasal_timestamps_REST) && ~isempty(sd.temporal_timestamps_REST)   % There are some sessions with zero spontaneous (platform stationary) saccades
        cfg_1.doPlot = 0;
        cfg_1.window = [-2 2];
        cfg_1.dt = 0.005;
        [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, S, sd.nasal_timestamps_REST);
        [mn2, edges] = histcounts(outputS_n, outputIT_n);
        plot(edges(1:end-1), mn2/cfg_1.dt/length(sd.nasal_timestamps_MOVING), 'LineWidth', cfg.LineWidth);
        amax = max(mn2/cfg_1.dt/length(sd.nasal_timestamps_MOVING));
        
        hold on
        [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, S, sd.temporal_timestamps_REST);
        [mt2, edges] = histcounts(outputS_t, outputIT_t);
        plot(edges(1:end-1), smoothdata(mt2/cfg_1.dt/length(sd.temporal_timestamps_MOVING)), 'LineWidth', cfg.LineWidth);
        bmax = max(smoothdata(mt2/cfg_1.dt/length(sd.temporal_timestamps_MOVING)));
        allmax = max([amax bmax]);
        c = axis;
        line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
        legend('nasal', 'temporal', '', 'FontSize', cfg.smallfont)
        set(gca, 'FontSize', cfg.FontSize)
        % set(gca, 'XTick', [])
        set(gca, 'XTick', [-2 2])
        axis([c(1) c(2) c(3) allmax])
    end
    text(NaN, NaN, strcat('#temp.=', sprintf('%0.0f', sd.numTemporal_stationary)), 'FontSize', 16, 'Units', 'normalized', 'Position', [.02 .95 0])
    text(NaN, NaN, strcat('#nasal=', sprintf('%0.0f', sd.numNasal_stationary)), 'FontSize', 16, 'Units', 'normalized', 'Position', [.02 .80 0])
    
end
%% #9 AHV PETH
subtightplot(6,6,9, [cfg.tightX cfg.tightY]); hold on
cfg_peth.window = [-2 2];
cfg_peth.mode = 'interp';
cfg_peth.dt = median(diff(sd.AHV.tvec));
out_nasal = TSDpeth(cfg_peth, sd.AHV, sd.nasal_timestamps_MOVING);
out_temporal = TSDpeth(cfg_peth, sd.AHV, sd.temporal_timestamps_MOVING);
plot(out_nasal, 'LineWidth', cfg.LineWidth); hold on
plot(out_temporal, 'LineWidth', cfg.LineWidth);
set(gca, 'FontSize', cfg.FontSize)
axis tight
c = axis;
axis([-2 2 c(3) c(4)])
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
legend('nasal', 'temporal', '', 'FontSize', cfg.smallfont)
set(gca, 'XTick', [-2 2])
title('AHV PETH', 'Units', 'normalized', 'Position', [0.5, 0.5, 0], 'FontSize', cfg.FontSize);

%% #10 HistISI
p = subtightplot(6,6,10, [cfg.tightX cfg.tightY]); hold on
[h, n] = HistISIsubplot(S.t{1});
HistISIsubplot(S.t{1});
c = axis;
[~, i] = max(h);
line([n(i) n(i)], [0 c(4)], 'color', 'k');
grid on
set(gca, 'TickDir', 'out', 'XLim', [cfg.histXmin cfg.histXmax], 'FontSize', cfg.FontSize)
title('HIST ISI', 'Units', 'normalized', 'Position', [0.5, 0.5, 0], 'FontSize', cfg.FontSize);
set(gca, 'XTick', [])

%% #10 laser-triggered EYE movement
subtightplot(6,6,11, [cfg.tightX cfg.tightY]); hold on
% out: tsd with PETH
% cfg options:
%
% if eye
%     if cfg.doLaser == 1
%         [~, ~, laser_on, laser_off, arraysize, ~] = SortBrainstemEventLabels3;
%         dur = mode(laser_off - laser_on);
%         
%         cfg_eye.window = [-2 2];
%         cfg_eye.dt = .1;
%         title('Laser-Aligned Eye Movement')
%         out = TSDpeth(cfg_eye, sd.tsdH, laser_on);   
%     end
% end
%% #11 BLANK

%% #12 BLANK

%% #13 MOVING SACCADE peth: narrow
if eye
    subtightplot(6,6,13, [cfg.tightX cfg.tightY]); hold on
    if ~isempty(sd.nasal_timestamps_MOVING) && ~isempty(sd.temporal_timestamps_MOVING)
        cfg_1.doPlot = 0;
        cfg_1.window = [-.2 .2];
        cfg_1.dt = 0.01;
        [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, S, sd.nasal_timestamps_MOVING);
        [mn3, edges] = histcounts(outputS_n, outputIT_n);
        plot(edges(1:end-1), mn3/cfg_1.dt/length(sd.nasal_timestamps_MOVING), 'LineWidth', cfg.LineWidth);
        amax = max( mn3/cfg_1.dt/length(sd.nasal_timestamps_MOVING));
        set(gca, 'FontSize', cfg.FontSize)
        
        hold on
        [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, S, sd.temporal_timestamps_MOVING);
        [mt3, edges] = histcounts(outputS_t, outputIT_t);
        plot(edges(1:end-1), mt3/cfg_1.dt/length(sd.temporal_timestamps_MOVING), 'LineWidth', cfg.LineWidth);
        bmax = max(mt3/cfg_1.dt/length(sd.temporal_timestamps_MOVING));
        allmax = max([amax bmax]);
        title('Sacc. peth MOVING')
        ylabel('FR (Hz)', 'FontSize', cfg.FontSize)
        set(gca, 'FontSize', cfg.FontSize)
        c = axis;
        line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
        axis([c(1) c(2) c(3) allmax])
    end
    
    %% #14 Stationary SACCADE peth: narrow
    subtightplot(6,6,14, [cfg.tightX cfg.tightY]); hold on
    
    if ~isempty(sd.nasal_timestamps_REST) && ~isempty(sd.temporal_timestamps_REST)
        cfg_1.doPlot = 0;
        cfg_1.window = [-.2 .2];
        cfg_1.dt = 0.01;
        [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, S, sd.nasal_timestamps_REST);
        [mn4, edges] = histcounts(outputS_n, outputIT_n);
        plot(edges(1:end-1), mn4/cfg_1.dt/length(sd.nasal_timestamps_MOVING), 'LineWidth', cfg.LineWidth);
        amax = max(mn4/cfg_1.dt/length(sd.nasal_timestamps_MOVING));
        hold on
        [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, S, sd.temporal_timestamps_REST);
        [mt4, edges] = histcounts(outputS_t, outputIT_t);
        plot(edges(1:end-1), mt4/cfg_1.dt/length(sd.temporal_timestamps_MOVING), 'LineWidth', cfg.LineWidth);
        bmax = max(mt4/cfg_1.dt/length(sd.temporal_timestamps_MOVING));
        allmax = max([amax bmax]);
        title('Sacc. STATIONARY')
        set(gca, 'FontSize', cfg.FontSize)
        xlabel('time peri saccade (s)')
        c = axis;
        line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
        axis([c(1) c(2) c(3) allmax])
    end
end
%% #15 Laser PETH
p = subtightplot(6,6,15, [cfg.tightX cfg.tightY]); hold on
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Get the correct timestamps for the different session events
if cfg.doLaser == 1
    [~, ~, laser_on, laser_off, arraysize, ~] = SortBrainstemEventLabels3;
    dur = mode(laser_off - laser_on);
    
    if ~isnan(laser_on)
        cfg_laser.window = [-1.1 2];
        cfg_laser.dt = 0.01;
        cfg_laser.binsize = cfg_laser.dt; % used for gaussian kernal.  select a small bin size for good time resolution
        cfg_laser.doPlot = 1;
        cfg_laser.doRaster = 1;
        cfg_laser.doBar = 0;
        [~, ~, ~, ~, ~] = SpikePETH_either(cfg_laser, sd.S, laser_on); hold on
        c = axis;
        rectangle(Position = [0, 0, dur, c(4)], FaceColor=[0 1 1], EdgeColor=[0 1 1])    % change the 3rd entry in rectangle to the diff of laser_off and laser_on
        [~, ~, ~, ~, ~] = SpikePETH_either(cfg_laser, sd.S, laser_on);  % repeating spikepeth here is a hack to get the background color 'in back'. otherwise it occludes the spikes.
        line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
        line([1 1], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
        set(gca, 'YTick', [])
        set(gca, 'FontSize', cfg.FontSize)
        title('Laser PETH', 'FontSize', cfg.FontSize)
        set(gca, 'XTick', [-1 0 1 2])
        hold on
        % add line for Ave Firing Rate, esp. when rasters are so dense that it turns black
        cfg_laser.window = [-1.1 2];
        cfg_laser.dt = 0.05;
        cfg_laser.binsize = cfg_laser.dt; % used for gaussian kernal.  select a small bin size for good time resolution
        cfg_laser.doPlot = 0;
        cfg_laser.doRaster = 0;
        cfg_laser.doBar = 0;
        [outputS_laser, ~, ~, outputIT_laser, ~] = SpikePETH_either(cfg_laser, sd.S, laser_on);
        m = histc(outputS_laser, outputIT_laser);
        yyaxis right
        if ~isempty(outputS_laser)
            plot(outputIT_laser(1:end-1),m(1:end-1)/cfg_laser.dt/length(laser_on), 'LineWidth', 4);   % JJS. 11/15/22. This is a hack. Last value of m is always zero (erroneously).
        end
    else
        disp('laser events not found')
        title('Laser PETH', 'FontSize', cfg.FontSize)
    end
end

%% #16 DUMMY Peth
p = subtightplot(6,6,16, [cfg.tightX cfg.tightY]); hold on
if cfg.doLaser == 1
    events_ts = LoadEvents([]); %#ok<*NASGU>
    wrapper = @(events_ts) strcmp(events_ts, 'ShutterSound On');
    A = cellfun(wrapper, sd.Events.label);
    sound_label = find(A); % index which label corresponds to 'ShutterSound On'
    if ~isempty(sound_label)
        sound_times = sd.Events.t{sound_label};
        sound_times = sound_times - sd.starttime;
        
        cfg_dummy.window = [-1.1 2];
        cfg_dummy.dt = 0.01;
        cfg_dummy.binsize = cfg_dummy.dt; % used for gaussian kernal.  select a small bin size for good time resolution
        cfg_dummy.doPlot = 1;
        cfg_dummy.doRaster = 1;
        cfg_dummy.doBar = 0;
        [~, ~, ~, ~, ~] = SpikePETH_either(cfg_dummy, sd.S, sound_times); hold on
        c = axis;
        rectangle(Position = [0, 0, dur, c(4)], FaceColor=[0 .9 .4], EdgeColor=[0 .9 .4])    % change the 3rd entry in rectangle to the diff of laser_off and laser_on
        [~, ~, ~, ~, ~] = SpikePETH_either(cfg_dummy, sd.S, sound_times);  % this is a hack to get the background color 'in back'. otherwise it occludes the spikes.
        line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
        line([1 1], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
        set(gca, 'YTick', [])
        set(gca, 'FontSize', cfg.FontSize)
        title('Dummy PETH', 'FontSize', cfg.FontSize)
        set(gca, 'XTick', [-1 0 1 2])
        set(gca, 'Ylabel', [])
        hold on
        
        % add line for Ave Firing Rate, esp. when rasters are so dense that it turns black
        cfg_dummy.dt = 0.05;
        cfg_dummy.doPlot = 0;
        cfg_dummy.doRaster = 0;
        cfg_dummy.doBar = 0;
        [outputS_dummy, ~, ~, outputIT_dummy, ~] = SpikePETH_either(cfg_dummy, sd.S, sound_times);
        m = histc(outputS_dummy, outputIT_dummy); %#ok<*HISTC>
        % m = histcounts(outputS, outputIT);
        yyaxis right
        plot(outputIT_dummy(1:end-1),m(1:end-1)/cfg_dummy.dt/length(sound_times), 'LineWidth', 4);   % JJS. 11/15/22. This is a hack. Last value of m is always zero (erroneously).
    else
        warning('this session lacks shutter sound timestamps.')
        title('Dummy PETH', 'FontSize', cfg.FontSize)
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
    end
else
    disp('Laser events, if they exist, skipped')
end

%% #17 Average Waveform
p = subtightplot(6,6,17, [cfg.tightX cfg.tightY]); hold on
% title('Average Waveform', 'FontSize', FontSize)
waveName = strcat(newID, '-wv.mat');
if exist(waveName)
    load(waveName)
    CO = get(gca,'ColorOrder');
    eIND = xrange(1:2:end,1);
    hold on
    %Extract AD2BitVoltConversioFactor from corresponding .ntt file
    TTpos = strfind(cellID, 'TT');
    TTnum = cellID(TTpos+2:TTpos+3);   % position of the numerals after 'TT'
    NT = FindFile(strcat('*', TTnum, '.ntt'));
    [~, b, c] = fileparts(NT);
    filename = strcat(b,c);
    
    H = ReadNewHeader(filename);
    L = find(~cellfun(@isempty,strfind(H, '-ADBitVolts')));
    [~,remain] = strtok(H{L});
    CONV = textscan(remain, '%f'); % conversion factor
    
    for iT = 1:4
        conv = CONV{1}(iT)*10^6;
        temp = xrange(:,iT);
        ind = temp(1:2:end);
        plot(xrange(:,iT), mWV(:,iT)*conv);
        h = errorbar(ind, mWV(eIND,iT)*conv, sWV(eIND,iT)*conv,'.');
        set(h, 'Color',CO(2*iT-1,:),'MarkerSize',1); % it uses every other color, for some reason
    end
    if exist(strcat(newID, '-ClusterQual.mat'))
        load(strcat(newID, '-ClusterQual.mat'))
        text(NaN, NaN, sprintf('Lratio %.2f \n', CluSep.Lratio), 'FontSize', 13, 'Units', 'normalized', 'Position', [.1 .85 0])
        text(NaN, NaN, sprintf('IsoD %.1f \n', CluSep.IsolationDist), 'FontSize', 13, 'Units', 'normalized', 'Position', [.6 .85 0])
    end
    
else
    disp('no wave file found. skipping average waveform')
end
set(gca, 'XTick', [])
set(gca, 'FontSize', cfg.FontSize)
axis tight
p.YAxisLocation = 'right';

%% #18 Session ID
subtightplot(6,6, 18, [cfg.tightX cfg.tightY]);
title(newID, 'Color', 'r', 'FontSize', 20)
set(gca, 'XTick', [])
set(gca, 'YTick', [])

%% #19:24 WHEEL SPEED and AHV
plot19 = subtightplot(6,6,19:24, [cfg.tightX cfg.tightY]);

if exist(WheelencoderCSC)
    yyaxis left
    plot(sd.speed.tvec, sd.speed.data, 'Color', 'm')
    set(gca, 'Xlim', [0 sd.speed.tvec(end)], 'FontSize', cfg.FontSize)
    ylabel('Wheel Speed (cm/s)', 'Color', 'm')
    yyaxis right
    plot(sd.AHV.tvec, sd.AHV.data)
    set(gca, 'XTick', [])
    ylabel('AHV')
    
    yyaxis left
    c = axis; hold on
    if cfg.doLaser
        if ~isnan(laser_on)
            for iLaser = 1:arraysize
                rectangle(Position=[laser_on(iLaser), 0, laser_off(iLaser) - laser_on(iLaser), c(4)], FaceColor=[0 1 1], EdgeColor=[0 1 1])
            end
        end
        if exist('sound_times')
            for iSound = 1:length(sound_times)
                rectangle(Position=[sound_times(iSound), 0, 1, c(4)], FaceColor=[0 .9 .4], EdgeColor=[0 .9 .4])  % JJS. 11/19/22. For now, shutter sound duration is 1 second.
            end
        end
    end
end

%%  #25:30 FIRING RATE and AHV
plot25 = subtightplot(6,6,25:30, [cfg.tightX cfg.tightY]); hold on
cfg_Q = []; cfg_Q.dt = 0.001; cfg_Q.gausswin_sd = 0.05;cfg_Q.smooth = 'gauss';
Q = MakeQfromS(cfg_Q, sd.S);
tvec = Q.tvec - Q.tvec(1);
yyaxis left
h = plot(Q.tvec, Q.data./cfg_Q.dt);
set(gca, 'Xlim', [0 tvec(end)], 'FontSize', cfg.FontSize)
ylabel('FR (Hz)', 'FontSize', cfg.FontSize)
yyaxis right
plot(sd.AHV.tvec, sd.AHV.data)
ylabel('AHV')
% set(gca, 'XTick', [])
yyaxis left
c = axis;
if cfg.doLaser == 1
    if ~isnan(laser_on)
        for iLaser = 1:arraysize
            rectangle(Position=[laser_on(iLaser), 0, laser_off(iLaser) - laser_on(iLaser), c(4)], FaceColor=[0 1 1], EdgeColor=[0 1 1])
        end
    end
    if exist('sound_times')
        for iSound = 1:length(sound_times)
            rectangle(Position=[sound_times(iSound), 0, 1, c(4)], FaceColor=[0 .9 .4], EdgeColor=[0 .9 .4])  % JJS. 11/19/22. For now, shutter sound duration is 1 second.
        end
    end
    g = plot(Q.tvec, Q.data./cfg_Q.dt, 'LineStyle', '-');
end
% g.FaceColor = h.FaceColor;

%% #31:36 HORIZONTAL EYE POSITIION and AHV

plot31 = subtightplot(6,6,31:36, [cfg.tightX cfg.tightY]);
if exist(strcat(sd.SSN, '-VT1_proc.mat'))
    load(strcat(sd.SSN, '-VT1_proc.mat'), 'pupil');         % load the output of facemap
    [~, b, c] = fileparts(FindFile('*VT1.smi'));
    fn = strcat(b,c);
    tvec_raw = read_smi(fn);
    tvec = tvec_raw - sd.starttime;
    if strcmp(sd.SSN, 'M281-2021-12-23')              % exception for this session where cheetah crashed and .smi is shorter than pupilH
        tvec = .02*(1:length( pupil{1}.com));
        tvec = tvec';
    end
    yyaxis left
    plot(tvec, pupil{1}.com(:,2), 'Color', 'k');         % pupil{1}.com(:,2)  is the horizontal eye position from facemap
    ylabel('Horiz. Eye Pos. (pixels)', 'FontSize', cfg.FontSize, 'Color', 'k')
    xlabel('Time (sec)', 'FontSize', cfg.FontSize)
    set(gca, 'FontSize', cfg.FontSize)
    
    yyaxis right
    plot(sd.AHV.tvec, sd.AHV.data)
    ylabel('AHV (deg./sec)')
    c = axis;
    axis([c(1) sd.AHV.tvec(end) c(3) c(4)]);
    %%
    linkaxes([plot19 plot25 plot31], 'x');
    
    yyaxis left
    c = axis; hold on
    if cfg.doLaser
        if ~isnan(laser_on)
            for iLaser = 1:arraysize
                rectangle(Position=[laser_on(iLaser), 0, laser_off(iLaser) - laser_on(iLaser), c(4)], FaceColor=[0 1 1], EdgeColor=[0 1 1])
            end
        end
        if exist('sound_times')
            for iSound = 1:length(sound_times)
                rectangle(Position=[sound_times(iSound), 0, 1, c(4)], FaceColor=[0 .9 .4], EdgeColor=[0 .9 .4])  % JJS. 11/19/22. For now, shutter sound duration is 1 second.
            end
        end
    end
else
    disp('no eyetracking data for this session')
end

toc
end

