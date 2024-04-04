function plot_AHVtuningMegaplot6(iCell, varargin)
% JJS.
% For plotting most/all of the relevant data for a single cell for a headfixed brainstem recording session.
% 2021-02-16. Added more elements, like platform orientation and eye position. Expanded from 3x6 subtightplot to 4x6.
% 2022-10-19. Added and modified some elements. Important changes are:
%       - addition of saccade peths
%       - addition of laser peth (if applicable)    [inlcude blue background coloring to indicate laser on period]
%       - collapsed FR plot into on fig.
%       - adopted subplottight instead of subtightplot to maximize space
%       - include wheel speed info into one of the subplots
%       [aspirational: add ability to include histology image in one of the subtightplot spaces]    imread.m
tic
clf
tightX = .025;
tightY = .02;
occthresh = 0.5;
smallfont = 8;
insetText = 18;
process_varargin(varargin)
% cd to data folder
SSN = HD_GetSSN; disp(SSN);
% miscellaney
clf
FontSize = 13;
histXmin = 0.01;
histXmax = 0.2;
LineWidth = 3;
subtractStartTime = 1;
% Load Spikes
cfg = [];
cfg.uint = '64';
spikefiles = FindFiles('*.t');
cfg.fc = {spikefiles{iCell}};
Sold = LoadSpikes(cfg);

if subtractStartTime == 1 % New cheetah versions have timestamps
    if exist('events_ts.mat') == 2
        load('events_ts.mat');
    else
        events_ts = LoadEvents([]);
    end
    wrapper = @(events_ts) strcmp(events_ts, 'Starting Recording');
    A = cellfun(wrapper, events_ts.label);
    Startindex = find(A); % index which label says 'Start Recording'
    starttime = events_ts.t{Startindex}(1); % use the very first start record time
    for iC = 1:length(Sold.t)
        if strcmp(SSN, 'M054-2020-03-03') == 1
            S.t{iC} = Sold.t{iC} - Sold.t{iC}(1);
        else
            S.t{iC} = Sold.t{iC} - starttime;  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
        end
    end
end
% get AHV Tuning Curve
cfg_AHV = [];
cfg_AHV.subsample_factor = 10;
[AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, S);
AHV_dt = median(diff(AHV_tsd.tvec));
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #1 plot scatterplot
p = subtightplot(6,6,1, [tightX tightY]);
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = AHV_dt;
cfg_Q.tvec_edges = AHV_tsd.tvec(1):AHV_dt:AHV_tsd.tvec(end);
F = MakeQfromS(cfg_Q, S); % convert to FR
F.data = F.data ./ cfg_Q.dt;
F_idx = nearest_idx3(AHV_tsd.tvec, F.tvec);
AHV_F = F.data(:,F_idx);
ymax = max(AHV_F);
set(gca, 'TickDir', 'out', 'FontSize', FontSize)
plot(AHV_tsd.data, AHV_F, '.', 'MarkerSize', .5); hold on
set(gca, 'Ylim', [0 ymax], 'FontSize', FontSize)

% Add Tuning Curve
plot(tc_out.binCenters, tc_out.tc, 'LineWidth', LineWidth, 'Color', 'k');
ylabel('FR (Hz)', 'FontSize', FontSize)
set(groot, 'DefaultLegendInterpreter', 'none')
title('AHV Tuning Curve')
h = get(gca, 'XLim');
% text(.75*h(1), 10, 'CW', 'FontSize', 12)
% text(.5*h(2), 10, 'CCW', 'FontSize', 12)
text(NaN, NaN, 'CW', 'FontSize', insetText, 'Units', 'normalized', 'Position', [.25 .85 0])
text(NaN, NaN, 'CCW', 'FontSize', insetText, 'Units', 'normalized', 'Position', [.6 .85 0])
p.XAxisLocation = 'top';
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #2 pupil TC
p = subtightplot(6,6,2, [tightX tightY]); hold on
try
    load(FindFile('*saccades-edited.mat'), 'tsdH') % tsdH is the pupil position variable
catch
    warning('cannot find saccade data')
end
% calculate Q matrix
tsdH_dt = median(diff(tsdH.tvec));
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = tsdH_dt;
cfg_Q.tvec_edges = tsdH.tvec(1):tsdH_dt:tsdH.tvec(end);
F = MakeQfromS(cfg_Q, S); % convert to FR
F.data = F.data ./ cfg_Q.dt;

% find FR corresponding to each AHV sample
F_idx = nearest_idx3(tsdH.tvec, F.tvec);
tsdH_F = F.data(:,F_idx);
plot(tsdH.data, tsdH_F(1,:), '.', 'MarkerSize', .5, 'color', [.8 .8 .8]); hold on

tsdH.data = tsdH.data'; % change the shape so that it is a "well-formed tsd" for tuning curves
cfg_tc = [];
cfg_tc.nBins = 50;
cfg_tc.binEdges = {linspace(-60, 60, 101)};
cfg_tc.occ_dt = median(diff(tsdH.tvec));
cfg_tc.minOcc = 10;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
tc_pupil = TuningCurves(cfg_tc, S, tsdH);
plot(tc_pupil.usr.binCenters(tc_pupil.occ_hist>occthresh), smoothdata(tc_pupil.tc(1,(tc_pupil.occ_hist>occthresh))), 'k', 'LineWidth', 3);
set(gca, 'FontSize', FontSize)
title('Pupil Position (pixels)')
% axis tight
c = axis;
axis([-45 45 c(3) c(4)])
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
text(-30, c(4)/2, 'nasal', 'FontSize', insetText, 'Units', 'normalized', 'Position', [.15 .85 0])
text(5, c(4)/2, 'temporal', 'FontSize', insetText, 'Units', 'normalized', 'Position', [.55 .85 0])
p.XAxisLocation = 'top';

%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #3 acf
p = subtightplot(6,6,3, [tightX tightY]);
cfg_acf = [];
cfg_acf.binsize = 0.001;
cfg_acf.max_t = 0.5;
cfg_acf.smooth = 1; % set to 1 to compute ccf on SDF, 0 for raw spikes
cfg_acf.gauss_w = 1; % width of Gaussian convolution window (in s)
cfg_acf.gauss_sd = 0.005; % SD of Gaussian convolution window (in s)
set(gca, 'TickDir', 'out', 'FontSize', FontSize)
[acf, tvec] = ccf(cfg_acf, S.t{1}, S.t{1});
midpoint = ceil(length(acf)./2);
acf(midpoint) = 0;
plot(tvec, acf, 'LineWidth', 1);
% xlabel('Time (sec)', 'FontSize', FontSize)
set(gca, 'xtick', [-.5 -.25 0 .25 .5], 'FontSize', FontSize); grid on;
title('Acorr')
p.XAxisLocation = 'top';
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #4 acf zoomed in
p = subtightplot(6,6,4, [tightX tightY]); hold on
cfg_acf = [];
cfg_acf.binsize = 0.001;
cfg_acf.max_t = 0.05;
cfg_acf.smooth = 1; % set to 1 to compute ccf on SDF, 0 for raw spikes
cfg_acf.gauss_w = 1; % width of Gaussian convolution window (in s)
cfg_acf.gauss_sd = 0.005; % SD of Gaussian convolution window (in s)
set(gca, 'TickDir', 'out', 'FontSize', FontSize)
[acf, tvec] = ccf(cfg_acf, S.t{1}, S.t{1});
midpoint = ceil(length(acf)./2);
acf(midpoint) = 0;
plot(tvec, acf, 'LineWidth', 1);
% xlabel('Time (sec)', 'FontSize', FontSize)
set(gca, 'xtick', [-.05 0 .05], 'FontSize', FontSize); grid on;
title('Acorr')
p.XAxisLocation = 'top';
%---------------------------------------------------------------------------------------------------------------------------------------------------
%% #4 Upper Right Hand Corner: histo image
p = subtightplot(6,6,[5:6 11:12 17:18], [tightX tightY]); hold on
[fc] = FindFiles(strcat(SSN, '*.t'));
[~, b, ~] = fileparts(fc);
if iscell(b)
    cellID = b{iCell};
else
    cellID = b;
end
otherID = cellID;
otherID(21) = '-';
[~, filenamejpg, ext] = fileparts(FindFiles(strcat(cellID, '*', 'histology', '.jpg')));
filenamejpg = strcat(filenamejpg, ext);
[~, filenamepng, ext] = fileparts(FindFiles(strcat(cellID, '*', 'histology', '.png')));
filenamepng = strcat(filenamepng, ext);
doShow = 1;
if ~isempty(filenamejpg)
    filename = filenamejpg;
elseif ~isempty(filenamepng)
    filename = filenamepng;
else
    doShow = 0;
    disp('no histology file found for this neuron')
end
if doShow == 1
    [X,~] = imread(filename);
    imshow(X)
end
%% #10 HistISI
p = subtightplot(6,6,10, [tightX tightY]); hold on
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
p = subtightplot(6,6,15, [tightX tightY]); hold on
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------
laser_on = events_ts.t{1}-events_ts.t{2};                % this needs to be fixed. laser_on and laser_off occur at variable positions 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------------------------------------------
cfg_laser.window = [-1.1 2];
cfg_laser.dt = 0.01;
cfg_laser.binsize = cfg_laser.dt; % used for gaussian kernal.  select a small bin size for good time resolution
cfg_laser.doPlot = 1;
cfg_laser.doRaster = 1;
cfg_laser.doBar = 0;
[~, ~, ~, ~, ~] = SpikePETH_either(cfg_laser, S, laser_on); hold on
c = axis;
rectangle(Position=[0,0,1,c(4)], FaceColor=[0 1 1], EdgeColor=[0 1 1])
[~, ~, ~, ~, ~] = SpikePETH_either(cfg_laser, S, laser_on);  % this is a hack to get the background color 'in back'. otherwise it occludes the spikes.
set(gca, 'YTick', [])
set(gca, 'FontSize', FontSize)
title('Laser PETH', 'FontSize', FontSize)
set(gca, 'XTick', [-1 0 1 2])

%% #16 Laser PETH
p = subtightplot(6,6,16, [tightX tightY]); hold on
title('Average Waveform', 'FontSize', FontSize)
waveName = strcat(otherID, '-wv.mat');
if exist(waveName)
    load(waveName)
    CO = get(gca,'ColorOrder');
    eIND = xrange(1:2:end,1);
    hold on
    for iT = 1:4
        temp = xrange(:,iT);
        ind = temp(1:2:end);
        plot(xrange(:,iT), mWV(:,iT));
        h = errorbar(ind, mWV(eIND,iT), sWV(eIND,iT),'.');
        set(h, 'Color',CO(2*iT-1,:),'MarkerSize',1); % it uses every other color, for some reason
    end
else
    disp('no wave file found. skipping average waveform')
end
set(gca, 'XTick', []);
set(gca, 'YTick', []);

% %% #17 Average neuron waveform (non-laser spikes)
% p = subtightplot(6,6,17, [tightX tightY]); hold on
% title('Average Waveform: non-laser', 'FontSize', FontSize)
% % p.XAxisLocation = 'top';
% set(gca, 'Xtick', [])
%
% %% #18 Average neuron waveform (laser-only spikes)
% p = subtightplot(6,6,18, [tightX tightY]); hold on
% title('Average Waveform: laser only', 'FontSize', FontSize)
% % p.XAxisLocation = 'top';
% set(gca, 'Xtick', [])

%% GET SACCADE INFO
[~, ~, ~, ~, nasal_timestamps_MOVING, temporal_timestamps_MOVING] = isolateManualSaccades();
[~, ~, nasal_timestamps_REST, temporal_timestamps_REST] = isolateStationarySaccades();

%% #7 MOVING SACCADE peth: wide
subtightplot(6,6,7, [tightX tightY]); hold on
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
legend('nasal', 'temporal', '', 'FontSize', smallfont)
% set(gca, 'XTick', [])
set(gca, 'XTick', [-2 2])
ylabel('FR (Hz)', 'FontSize', FontSize)

%% #8 Stationary SACCADE peth: wide
subtightplot(6,6,8, [tightX tightY]); hold on

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
legend('nasal', 'temporal', '', 'FontSize', smallfont)
set(gca, 'FontSize', FontSize)
% set(gca, 'XTick', [])
set(gca, 'XTick', [-2 2])
axis([c(1) c(2) c(3) allmax])

%% #13 MOVING SACCADE peth: narrow
subtightplot(6,6,13, [tightX tightY]); hold on
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

%% #14 Stationary SACCADE peth: narrow
subtightplot(6,6,14, [tightX tightY]); hold on
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

%% #9 AHV PETH
subtightplot(6,6,9, [tightX tightY]); hold on
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
legend('nasal', 'temporal', '', 'FontSize', smallfont)
set(gca, 'XTick', [-2 2])
t = title('AHV PETH', 'Units', 'normalized', 'Position', [0.5, 0.5, 0], 'FontSize', FontSize);


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
plot7 = subtightplot(6,6,25:30, [tightX tightY]); hold on
cfg_Q = []; cfg_Q.dt = 0.001; cfg_Q.gausswin_sd = 0.05;cfg_Q.smooth = 'gauss';
Q = MakeQfromS(cfg_Q, S);
tvec = Q.tvec - Q.tvec(1);
yyaxis left
plot(Q.tvec, Q.data./cfg_Q.dt)
set(gca, 'Xlim', [0 tvec(end)], 'FontSize', FontSize)
ylabel('FR (Hz)', 'FontSize', FontSize)
yyaxis right
plot(AHV_tsd.tvec, AHV_tsd.data)
set(gca, 'XTick', [])

%% #19:24 WHEEL SPEED and AHV
plot8 = subtightplot(6,6,19:24, [tightX tightY]);

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

%% #31:36 HORIZONTAL EYE POSITIION and AHV
plot9 = subtightplot(6,6,31:36, [tightX tightY]);
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
else
end
toc


% title(Sold.label{iCell}, 'Color', 'r')
% t = title('this is my title', 'Units', 'normalized', 'Position', [0.5, 0.75, 0]);
% t.Color = 'r'; t.FontSize = 10; % with this you can change color, font name and size