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
cfg_def.doLaser = 0;        % Choose 0 to exclude laser events. Choose 1 to inlude laser events. This is a temporary measure b/c there are bugs in some sessions for laser event time extraction
cfg_def.FontSize = 13;      % General font size
cfg_def.histXmin = 0.01;    % Axis size for the HISTISI plot
cfg_def.histXmax = 0.2;
cfg_def.LineWidth = 3;      % General line width for plotting 

cfg = ProcessConfig2(cfg_def, cfg_in);

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
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #1 Firing Rate x AHV (with Tuning Curve) 
p = subtightplot(6,6,1, [cfg.tightX cfg.tightY]);
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = sd.AHV_dt;
cfg_Q.tvec_edges = sd.AHV.tvec(1): sd.AHV_dt: sd.AHV.tvec(end);
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
plot(tc_out.usr.binCenters, smoothdata(tc_out.tc(iCell,:)), 'LineWidth', cfg.LineWidth, 'Color', 'k');
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
p = subtightplot(6,6,2, [cfg.tightX cfg.tightY]); hold on
    % calculate Q matrix
    cfg_Q = [];
    cfg_Q.smooth = 'gauss';
    cfg_Q.gausswin_sd = 0.05;
    cfg_Q.dt = sd.tsdH_dt;
    cfg_Q.tvec_edges = sd.tsdH.tvec(1): sd.tsdH_dt: sd.tsdH.tvec(end);
    F = MakeQfromS(cfg_Q, sd.S); % convert to FR
    F.data = F.data(iCell,:) ./ cfg_Q.dt;
    
    % find FR corresponding to each AHV sample
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
    plot(tc_pupil.usr.binCenters(tc_pupil.occ_hist > cfg.occthresh), smoothdata(tc_pupil.tc(1,(tc_pupil.occ_hist > cfg.occthresh))), 'k', 'LineWidth', 3);
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












%% #10 HistISI
p = subtightplot(6,6,10, [cfg.tightX cfg.tightY]); hold on
[h, n] = HistISIsubplot(S.t{1});
HistISIsubplot(S.t{1});
c = axis;
[~, i] = max(h);
line([n(i) n(i)], [0 c(4)], 'color', 'k');
grid on
set(gca, 'TickDir', 'out', 'XLim', [histXmin histXmax], 'FontSize', FontSize)
% xlabel('Time (sec)', 'FontSize', FontSize)
% title('HistISI')
t = title('HIST ISI', 'Units', 'normalized', 'Position', [0.5, 0.5, 0], 'FontSize', FontSize);
% p.XAxisLocation = 'top';
set(gca, 'XTick', [])

%% #15 Laser PETH
p = subtightplot(6,6,15, [cfg.tightX cfg.tightY]); hold on
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Get the correct timestamps for the different session events
if cfg.doLaser == 1
    [start_time, stop_time, laser_on, laser_off, bit0, bit4, arraysize] = SortBrainstemEventLabels;
    % [start_time, stop_time, laser_on, laser_off, bit0, bit4] = SortBrainstemEventLabels2;
    
    if ~isnan(laser_on)
        %     dur = mode(laser_off(1:arraysize) - laser_on(1:arraysize));
        dur = mode(laser_off - laser_on);
        
        cfg_laser.window = [-1.1 2];
        cfg_laser.dt = 0.01;
        cfg_laser.binsize = cfg_laser.dt; % used for gaussian kernal.  select a small bin size for good time resolution
        cfg_laser.doPlot = 1;
        cfg_laser.doRaster = 1;
        cfg_laser.doBar = 0;
        [~, ~, ~, ~, ~] = SpikePETH_either(cfg_laser, S, laser_on); hold on
        c = axis;
        rectangle(Position = [0, 0, dur, c(4)], FaceColor=[0 1 1], EdgeColor=[0 1 1])    % change the 3rd entry in rectangle to the diff of laser_off and laser_on
        [outputS, ~, ~, outputIT, cfg] = SpikePETH_either(cfg_laser, S, laser_on);  % this is a hack to get the background color 'in back'. otherwise it occludes the spikes.
        line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
        line([1 1], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
        set(gca, 'YTick', [])
        set(gca, 'FontSize', FontSize)
        title('Laser PETH', 'FontSize', FontSize)
        set(gca, 'XTick', [-1 0 1 2])
        hold on
        % add line for Ave Firing Rate, esp. when rasters are so dense that it turns black
        cfg_laser.window = [-1.1 2];
        cfg_laser.dt = 0.05;
        cfg_laser.binsize = cfg_laser.dt; % used for gaussian kernal.  select a small bin size for good time resolution
        cfg_laser.doPlot = 0;
        cfg_laser.doRaster = 0;
        cfg_laser.doBar = 0;
        [outputS_laser, ~, ~, outputIT_laser, cfg_out_laser] = SpikePETH_either(cfg_laser, S, laser_on);
        m = histc(outputS_laser, outputIT_laser);
        yyaxis right
        plot(outputIT_laser(1:end-1),m(1:end-1)/cfg_laser.dt/length(laser_on), 'LineWidth', 4);   % JJS. 11/15/22. This is a hack. Last value of m is always zero (erroneously).
    else
        disp('laser events not found')
        title('Laser PETH', 'FontSize', FontSize)
    end
    %% #16 DUMMY Peth
    p = subtightplot(6,6,16, [cfg.tightX cfg.tightY]); hold on
    
    wrapper = @(events_ts) strcmp(events_ts, 'ShutterSound On');
    A = cellfun(wrapper, events_ts.label);
    sound_label = find(A); % index which label corresponds to 'ShutterSound On'
    if ~isempty(sound_label)
        sound_times = events_ts.t{sound_label};
        sound_times = sound_times - start_time;
        
        cfg_dummy.window = [-1.1 2];
        cfg_dummy.dt = 0.01;
        cfg_dummy.binsize = cfg_dummy.dt; % used for gaussian kernal.  select a small bin size for good time resolution
        cfg_dummy.doPlot = 1;
        cfg_dummy.doRaster = 1;
        cfg_dummy.doBar = 0;
        [~, ~, ~, ~, ~] = SpikePETH_either(cfg_dummy, S, sound_times); hold on
        c = axis;
        rectangle(Position = [0, 0, dur, c(4)], FaceColor=[0 .9 .4], EdgeColor=[0 .9 .4])    % change the 3rd entry in rectangle to the diff of laser_off and laser_on
        [outputS_dummy, ~, ~, outputIT_dummy, cfg_out_dummy] = SpikePETH_either(cfg_dummy, S, sound_times);  % this is a hack to get the background color 'in back'. otherwise it occludes the spikes.
        line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
        line([1 1], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'w')
        set(gca, 'YTick', [])
        set(gca, 'FontSize', FontSize)
        title('Dummy PETH', 'FontSize', FontSize)
        set(gca, 'XTick', [-1 0 1 2])
        set(gca, 'Ylabel', [])
        hold on
        
        % add line for Ave Firing Rate, esp. when rasters are so dense that it turns black
        cfg_dummy.dt = 0.05;
        cfg_dummy.doPlot = 0;
        cfg_dummy.doRaster = 0;
        cfg_dummy.doBar = 0;
        [outputS_dummy, ~, ~, outputIT_dummy, cfg_out_dummy] = SpikePETH_either(cfg_dummy, S, sound_times);
        m = histc(outputS_dummy, outputIT_dummy);
        % m = histcounts(outputS, outputIT);
        yyaxis right
        plot(outputIT_dummy(1:end-1),m(1:end-1)/cfg_dummy.dt/length(sound_times), 'LineWidth', 4);   % JJS. 11/15/22. This is a hack. Last value of m is always zero (erroneously).
        
    else
        warning('this session lacks shutter sound timestamps.')
        title('Dummy PETH', 'FontSize', FontSize)
    end
else
    disp('Laser events, if they exist, skipped')
end

%% #16 Average Waveform
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
set(gca, 'FontSize', FontSize)
axis tight
p.YAxisLocation = 'right';


%% GET SACCADE INFO
[~, ~, ~, ~, nasal_timestamps_MOVING, temporal_timestamps_MOVING] = isolateManualSaccades();
[~, ~, nasal_timestamps_REST, temporal_timestamps_REST] = isolateStationarySaccades();

numNasal_moving = length(nasal_timestamps_MOVING); if numNasal_moving == 1; numNasal_moving = 0; end  % if this array is a single NaN, then there are in fact no saccades.
numNasal_stationary = length(nasal_timestamps_REST); if numNasal_stationary == 1; numNasal_stationary = 0; end

numTemporal_moving = length(temporal_timestamps_MOVING); if numTemporal_moving == 1; numTemporal_moving = 0; end
numTemporal_stationary = length(temporal_timestamps_REST); if numTemporal_stationary == 1; numTemporal_stationary = 0; end

%% #7 MOVING SACCADE peth: wide
subtightplot(6,6,7, [cfg.tightX cfg.tightY]); hold on
% t = title('Evoked Saccade PETH', 'FontSize', FontSize, 'Units', 'normalized', 'Position', [.5 .1 0]);
cfg_1.doPlot = 0;
cfg_1.window = [-2 2];
cfg_1.dt = 0.05;
[outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, S, nasal_timestamps_MOVING);
[mn1, edges] = histcounts(outputS_n, outputIT_n);
plot(edges(1:end-1), mn1/cfg_1.dt/length(nasal_timestamps_MOVING), 'LineWidth', LineWidth); % this is a hack. should replace binedgges with bincenters

hold on
[outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, S, temporal_timestamps_MOVING);
[mt1, edges] = histcounts(outputS_t, outputIT_t);
plot(edges(1:end-1), mt1/cfg_1.dt/length(temporal_timestamps_MOVING), 'LineWidth', LineWidth);
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
set(gca, 'FontSize', FontSize)
legend('nasal', 'temporal', '', 'FontSize', cfg.smallfont)
% set(gca, 'XTick', [])
set(gca, 'XTick', [-2 2])
ylabel('FR (Hz)', 'FontSize', FontSize)

% text(NaN, NaN, sprintf('Lratio %.2f \n', CluSep.Lratio), 'FontSize', 10, 'Units', 'normalized', 'Position', [.2 .85 0])
text(NaN, NaN, strcat('#temp.=', sprintf('%0.0f', numTemporal_moving)), 'FontSize', 16, 'Units', 'normalized', 'Position', [.02 .95 0])
text(NaN, NaN, strcat('#nasal=', sprintf('%0.0f', numNasal_moving)), 'FontSize', 16, 'Units', 'normalized', 'Position', [.02 .80 0])


%% #8 Stationary SACCADE peth: wide
subtightplot(6,6,8, [cfg.tightX cfg.tightY]); hold on
if ~isempty(nasal_timestamps_REST) && ~isempty(temporal_timestamps_REST)
    cfg_1.doPlot = 0;
    cfg_1.window = [-2 2];
    cfg_1.dt = 0.005;
    [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, S, nasal_timestamps_REST);
    [mn2, edges] = histcounts(outputS_n, outputIT_n);
    plot(edges(1:end-1), mn2/cfg_1.dt/length(nasal_timestamps_MOVING), 'LineWidth', LineWidth);
    amax = max(mn2/cfg_1.dt/length(nasal_timestamps_MOVING));
    
    hold on
    [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, S, temporal_timestamps_REST);
    [mt2, edges] = histcounts(outputS_t, outputIT_t);
    plot(edges(1:end-1), smoothdata(mt2/cfg_1.dt/length(temporal_timestamps_MOVING)), 'LineWidth', LineWidth);
    bmax = max(smoothdata(mt2/cfg_1.dt/length(temporal_timestamps_MOVING)));
    allmax = max([amax bmax]);
    c = axis;
    line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
    legend('nasal', 'temporal', '', 'FontSize', cfg.smallfont)
    set(gca, 'FontSize', FontSize)
    % set(gca, 'XTick', [])
    set(gca, 'XTick', [-2 2])
    axis([c(1) c(2) c(3) allmax])
end
text(NaN, NaN, strcat('#temp.=', sprintf('%0.0f', numTemporal_stationary)), 'FontSize', 16, 'Units', 'normalized', 'Position', [.02 .95 0])
text(NaN, NaN, strcat('#nasal=', sprintf('%0.0f', numNasal_stationary)), 'FontSize', 16, 'Units', 'normalized', 'Position', [.02 .80 0])

%% #13 MOVING SACCADE peth: narrow
subtightplot(6,6,13, [cfg.tightX cfg.tightY]); hold on

if ~isempty(nasal_timestamps_MOVING) && ~isempty(temporal_timestamps_MOVING)
    cfg_1.doPlot = 0;
    cfg_1.window = [-.2 .2];
    cfg_1.dt = 0.01;
    [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, S, nasal_timestamps_MOVING);
    [mn3, edges] = histcounts(outputS_n, outputIT_n);
    plot(edges(1:end-1), mn3/cfg_1.dt/length(nasal_timestamps_MOVING), 'LineWidth', LineWidth);
    amax = max( mn3/cfg_1.dt/length(nasal_timestamps_MOVING));
    set(gca, 'FontSize', FontSize)
    
    hold on
    [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, S, temporal_timestamps_MOVING);
    [mt3, edges] = histcounts(outputS_t, outputIT_t);
    plot(edges(1:end-1), mt3/cfg_1.dt/length(temporal_timestamps_MOVING), 'LineWidth', LineWidth);
    bmax = max(mt3/cfg_1.dt/length(temporal_timestamps_MOVING));
    allmax = max([amax bmax]);
    title('Sacc. peth MOVING')
    ylabel('FR (Hz)', 'FontSize', FontSize)
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
    axis([c(1) c(2) c(3) allmax])
end

%% #14 Stationary SACCADE peth: narrow
subtightplot(6,6,14, [cfg.tightX cfg.tightY]); hold on

if ~isempty(nasal_timestamps_REST) && ~isempty(temporal_timestamps_REST)
    cfg_1.doPlot = 0;
    cfg_1.window = [-.2 .2];
    cfg_1.dt = 0.01;
    [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, S, nasal_timestamps_REST);
    [mn4, edges] = histcounts(outputS_n, outputIT_n);
    plot(edges(1:end-1), mn4/cfg_1.dt/length(nasal_timestamps_MOVING), 'LineWidth', LineWidth);
    amax = max(mn4/cfg_1.dt/length(nasal_timestamps_MOVING));
    hold on
    [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, S, temporal_timestamps_REST);
    [mt4, edges] = histcounts(outputS_t, outputIT_t);
    plot(edges(1:end-1), mt4/cfg_1.dt/length(temporal_timestamps_MOVING), 'LineWidth', LineWidth);
    bmax = max(mt4/cfg_1.dt/length(temporal_timestamps_MOVING));
    allmax = max([amax bmax]);
    title('Sacc. STATIONARY')
    set(gca, 'FontSize', FontSize)
    xlabel('time peri saccade (s)')
    c = axis;
    line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
    axis([c(1) c(2) c(3) allmax])
end

%% #9 AHV PETH
subtightplot(6,6,9, [cfg.tightX cfg.tightY]); hold on
cfg_peth.window = [-2 2];
cfg_peth.mode = 'interp';
cfg_peth.dt = median(diff(AHV_tsd.tvec));
out_nasal = TSDpeth(cfg_peth, AHV_tsd, nasal_timestamps_MOVING);
out_temporal = TSDpeth(cfg_peth, AHV_tsd, temporal_timestamps_MOVING);
plot(out_nasal, 'LineWidth', LineWidth); hold on
plot(out_temporal, 'LineWidth', LineWidth);
%                 ylabel('AHV (deg/s)')
set(gca, 'FontSize', FontSize)
axis tight
c = axis;
axis([-2 2 c(3) c(4)])
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
legend('nasal', 'temporal', '', 'FontSize', cfg.smallfont)
set(gca, 'XTick', [-2 2])
t = title('AHV PETH', 'Units', 'normalized', 'Position', [0.5, 0.5, 0], 'FontSize', FontSize);
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #10 laser-triggered EYE movement
subtightplot(6,6,11, [cfg.tightX cfg.tightY]); hold on
% out: tsd with PETH
% cfg options:
%
% cfg_def.window = [-2 2]; % start and end times of window (in s)
% cfg_def.dt = []; % time step, specify this for 'interp' mode
% cfg_def.mode = 'raw'; % 'raw' or 'interp'; if 'interp', need to specify cfg_def.dt
% cfg_def.interp_mode = 'linear';

cfg_eye.window = [-2 2];
cfg_eye.dt = .1;
% out = TSDpeth(cfg_eye, tsdH, laser_on);



% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%  #25:30 FIRING RATE and AHV
plot7 = subtightplot(6,6,25:30, [cfg.tightX cfg.tightY]); hold on
cfg_Q = []; cfg_Q.dt = 0.001; cfg_Q.gausswin_sd = 0.05;cfg_Q.smooth = 'gauss';
Q = MakeQfromS(cfg_Q, S);
tvec = Q.tvec - Q.tvec(1);
yyaxis left
h = plot(Q.tvec, Q.data./cfg_Q.dt);
set(gca, 'Xlim', [0 tvec(end)], 'FontSize', FontSize)
ylabel('FR (Hz)', 'FontSize', FontSize)
yyaxis right
plot(AHV_tsd.tvec, AHV_tsd.data)
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

%% #19:24 WHEEL SPEED and AHV
plot8 = subtightplot(6,6,19:24, [cfg.tightX cfg.tightY]);

encoderCSC = strcat(SSN, '-CSC34.Ncs');
if exist(encoderCSC)
    yyaxis left
    updownTSD = getQEupdown([]);
    state_tsd = ConvertQEUpDownToState(updownTSD);
    [~, wheel_tsd] = ConvertQEStatesToAngle([], state_tsd);
    [d, speed, cfg_w] = ConvertWheeltoSpeed([], wheel_tsd);
    speed.tvec = speed.tvec - starttime;
    plot(speed.tvec, -speed.data)    % IMPORTANT: there is a sign change here b/c the voltage change in response to the mouse moving forward is negative. When the mouse moves 'forward', the wheel moves 'backward'.
    set(gca, 'Xlim', [0 tvec(end)], 'FontSize', FontSize)
    ylabel('Wheel Speed (cm/s)')
    yyaxis right
    plot(AHV_tsd.tvec, AHV_tsd.data)
    set(gca, 'XTick', [])
    
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
    plot(speed.tvec, -speed.data)
    
    
    %% #17   Wheel speed Tuning Curve
    p = subtightplot(6,6,5, [cfg.tightX cfg.tightY]); hold on
    % calculate Q matrix
    
    speed_dt = median(diff(tsdH.tvec));  % tsdH is the pupil position variable
    cfg_Q = [];
    cfg_Q.smooth = 'gauss';
    cfg_Q.gausswin_sd = 0.05;
    cfg_Q.dt = tsdH_dt;
    cfg_Q.tvec_edges = speed.tvec(1): speed_dt: speed.tvec(end);
    F = MakeQfromS(cfg_Q, S); % convert to FR
    F.data = F.data ./ cfg_Q.dt;
    
    speed.data = -speed.data; % we want forward motion to be displayed as a positive velocity
    % find FR corresponding to each AHV sample
    F_idx = nearest_idx3(speed.tvec, F.tvec);
    tsdH_F = F.data(:,F_idx);
    yyaxis right
    
    z = speed.data > cfg.speedthresh | speed.data < -cfg.speedthresh;
    plot(speed.data(z), tsdH_F(1,z), '.', 'MarkerSize', .5, 'color', [.8 .8 .8]); hold on
    %     if ~isempty(h) ==1
    h = lsline;
    set(h(1), 'Color', 'k')
    set(h(1), 'LineWidth', 2)
    %     end
    %     Fit = polyfit(h(1).XData, h(1).YData, 1);
    speeddata = speed.data(z)';
    speeddata(:,2) = ones(length(speeddata), 1);
    [b,bint,r,rint,stats] = regress(tsdH_F(1,z)', speeddata);
    c = axis;
    line([-cfg.speedthresh -cfg.speedthresh], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
    line([cfg.speedthresh cfg.speedthresh], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
    text(NaN, NaN, strcat('Rsq =', sprintf('%0.2f', stats(1))), 'FontSize', 12, 'Units', 'normalized', 'Position', [.55 .85 0])
    
    cfg_tc = [];
    cfg_tc.nBins = 50;
    cfg_tc.binEdges = {linspace(-5, 30, 101)};
    cfg_tc.occ_dt = median(diff(speed.tvec));
    cfg_tc.minOcc = 1;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
    tc_speed = TuningCurves(cfg_tc, S, speed);
    plot(tc_speed.usr.binCenters(tc_speed.occ_hist>cfg.occthresh), smoothdata(tc_speed.tc(1,(tc_speed.occ_hist>cfg.occthresh))), 'k', 'LineWidth', 3);
    axis tight
    ylabel('FR (Hz)', 'FontSize', FontSize)
    title('Wheel Speed Tuning Curve', 'FontSize', FontSize)
    yyaxis left
    set(gca, 'YTick', [])
    yyaxis right
    grid on
    set(gca, 'FontSize', FontSize)
    p.XAxisLocation = 'top';
    axis tight
    
    %% #31:36 WHEEL speed vs. AHV scatterplot
    plot6 = subtightplot(6,6,6, [cfg.tightX cfg.tightY]);
    
    X = AHV_tsd.data > cfg.ahv_thresh | AHV_tsd.data < -cfg.ahv_thresh;
    % Y = speed.data > cfg.speedthresh | speed.data < -cfg.speedthresh;
    
    Z_idx = nearest_idx3(AHV_tsd.tvec(X), speed.tvec);   % this did not work: Z_idx = nearest_idx3(AHV_tsd.tvec(X), speed.tvec(Y));
    SPEED = tsd(speed.tvec(Z_idx,:), speed.data(:,Z_idx));
    plot(AHV_tsd.data(X), SPEED.data, '.', 'MarkerSize', 1);
    axis tight
    xlabel('AHV', 'FontSize', FontSize)
    ylabel('Wheel Speed', 'FontSize', FontSize)
    plot6.YAxisLocation = 'right';
    plot6.XAxisLocation = 'top';
    set(gca, 'FontSize', FontSize)
    h = lsline;
    set(h(1), 'Color', 'k')
    set(h(1), 'LineWidth', 2)
    c = axis;
    line([c(1) c(2)], [0 0], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
    line([-cfg.ahv_thresh -cfg.ahv_thresh], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
    line([cfg.ahv_thresh cfg.ahv_thresh], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
    
    SPEEDregress = SPEED.data';
    SPEEDregress(:,2) = ones(length(SPEED.data), 1);
    [b,bint,r,rint,stats] = regress(AHV_tsd.data(X)', SPEEDregress);
    text(NaN, NaN, strcat('Rsq =', sprintf('%0.2f', stats(1))), 'FontSize', 12, 'Units', 'normalized', 'Position', [.05 .85 0])
    
else
    disp('WHEEL Encoder CSC not found')
    p = subtightplot(6,6,5, [cfg.tightX cfg.tightY]); hold on
    title('Wheel Speed Tuning Curve', 'FontSize', FontSize)
    plot6 = subtightplot(6,6,6, [cfg.tightX cfg.tightY]);
    
    xlabel('AHV', 'FontSize', FontSize)
    ylabel('Wheel Speed', 'FontSize', FontSize)
    plot6.YAxisLocation = 'right';
    plot6.XAxisLocation = 'top';
end

%% #31:36 HORIZONTAL EYE POSITIION and AHV
plot9 = subtightplot(6,6,31:36, [cfg.tightX cfg.tightY]);
SSN = HD_GetSSN;
if exist(strcat(SSN, '-VT1_proc.mat'))
    load(strcat(SSN, '-VT1_proc.mat'), 'pupil');         % load the output of facemap
    [a, b, c] = fileparts(FindFile('*VT1.smi'));
    fn = strcat(b,c);
    tvec_raw = read_smi(fn);
    tvec = tvec_raw - starttime;
    if strcmp(SSN, 'M281-2021-12-23')              % exception for this session where cheetah crashed and .smi is shorter than pupilH
        tvec = .02*(1:length( pupil{1}.com));
        tvec = tvec';
    end
    yyaxis left
    plot(tvec, pupil{1}.com(:,2));         % pupil{1}.com(:,2)  is the horizontal eye position from facemap
    ylabel('Horiz. Eye Pos. (pixels)', 'FontSize', FontSize)
    xlabel('Time (sec)', 'FontSize', FontSize)
    set(gca, 'FontSize', FontSize)
    
    yyaxis right
    plot(AHV_tsd.tvec, AHV_tsd.data)
    ylabel('AHV (deg./sec)')
    c = axis;
    axis([c(1) AHV_tsd.tvec(end) c(3) c(4)]);
    %%
    linkaxes([plot7 plot8 plot9], 'x');
    
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
    plot(tvec, pupil{1}.com(:,2), 'LineStyle', '-');
else
    disp('no eyetracking data for this session')
end

%% #18 Session ID
p = subtightplot(6,6, 18, [cfg.tightX cfg.tightY]);
title(newID, 'Color', 'r', 'FontSize', 20)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
toc



end

