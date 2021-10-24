function plot_AHVtuningMegaplot3(iCell, varargin)
% JJS. 2020.
% For plotting most/all of the relevant data for a single cell for a headfixed brainstem recording session.
% 2021-02-16. Added more elements, like platform orientation and eye position. Expanded from 3x6 subplot to 4x6.

process_varargin(varargin)
% cd to data folder
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
% cfg.fc = {spikefiles}; 
cfg.fc = {spikefiles{iCell}};
S = LoadSpikes(cfg);
% S = LoadSpikesJeff;

if subtractStartTime == 1 % New cheetah versions have timestamps
    events_ts = LoadEvents([]);
    wrapper = @(events_ts) strcmp(events_ts, 'Starting Recording');
    A = cellfun(wrapper, events_ts.label); 
    Startindex = find(A); % index which label says 'Start Recording' 
    starttime = events_ts.t{Startindex}(1); % use the very first start record time
    wrapper2 = @(events_ts) strcmp(events_ts, 'Stopping Recording');
    B = cellfun(wrapper2, events_ts.label);
    EndIndex = find(B); 
    endtime = events_ts.t{EndIndex};
%     assert(strcmp(events_ts.label{1}, 'Starting Recording')==1)
    for iC = 1:length(S.t)
        S.t{iC} = S.t{iC} - starttime;  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
    end
end
% get AHV Tuning Curve
cfg_AHV = [];
cfg_AHV.subsample_factor = 10;
[AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, S);
AHV_dt = median(diff(AHV_tsd.tvec));

%% #2 plot scatterplot
subplot(4,6,2)
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = AHV_dt;
cfg_Q.tvec_edges = AHV_tsd.tvec(1):AHV_dt:AHV_tsd.tvec(end);
F = MakeQfromS(cfg_Q, S); % convert to FR
% convert to FR
F.data = F.data ./ cfg_Q.dt;
% find FR corresponding to each AHV sample
F_idx = nearest_idx3(AHV_tsd.tvec, F.tvec);
AHV_F = F.data(:,F_idx);
ymax = max(AHV_F);

set(gca, 'TickDir', 'out', 'FontSize', FontSize)
plot(AHV_tsd.data, AHV_F, '.', 'MarkerSize', .5);
set(gca, 'Ylim', [0 ymax], 'FontSize', FontSize)
xlabel('AHV (deg./sec)', 'FontSize', FontSize)
ylabel('FR (Hz)', 'FontSize', FontSize)
title('Scatterplot')
h = get(gca, 'XLim');
text(.75*h(1), 10, 'CW', 'FontSize', 12)
text(.5*h(2), 10, 'CCW', 'FontSize', 12)

%% #1 plot tuning curves
subplot(4,6,1)
fc = FindFiles('*.t', 'CheckSubdirs', 0);
[a, b, c] = fileparts(fc{iCell});
plot(tc_out.usr.binCenters, tc_out.tc, 'LineWidth', LineWidth);
set(gca, 'TickDir', 'out', 'FontSize', FontSize)
set(gca, 'Ylim', [0 ymax])
xlabel('AHV (deg./sec)', 'FontSize', FontSize)
ylabel('FR (Hz)', 'FontSize', FontSize)
set(groot, 'DefaultLegendInterpreter', 'none')
title('Tuning Curve')
h = get(gca, 'XLim');
text(.75*h(1), 10, 'CW', 'FontSize', 12)
text(.5*h(2), 10, 'CCW', 'FontSize', 12)


%% #3 acf
subplot(4,6,3)
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
xlabel('Time (sec)', 'FontSize', FontSize)
set(gca, 'xtick', [-.5 -.25 0 .25 .5], 'FontSize', FontSize); grid on;
title('Acorr')


%% #4 acf zoomed in
subplot(4,6,4); hold on
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
xlabel('Time (sec)', 'FontSize', FontSize)
set(gca, 'xtick', [-.05 0 .05], 'FontSize', FontSize); grid on;
title('Acorr')


%% #5 HistISI
subplot(4,6,5); hold on
[h, n] = HistISIsubplot(S.t{1});
HistISIsubplot(S.t{1});
c = axis;
[~, i] = max(h);
line([n(i) n(i)], [0 c(4)], 'color', 'k');
grid on
set(gca, 'TickDir', 'out', 'XLim', [histXmin histXmax], 'FontSize', FontSize)
xlabel('Time (sec)', 'FontSize', FontSize)
title('HistISI')
% ylabel('Count')


%% #6 tbd
subplot(4,6,6); hold on






%%  #7 Firing Rate
plot7 = subplot(4,6,7:12); hold on
cfg_Q = []; cfg_Q.dt = 0.001; cfg_Q.gausswin_sd = 0.05;cfg_Q.smooth = 'gauss';
Q = MakeQfromS(cfg_Q, S);
tvec = Q.tvec - Q.tvec(1);
yyaxis left
plot(tvec, Q.data./cfg_Q.dt)
maxFR = max(Q.data./cfg_Q.dt);
set(gca, 'Xlim', [0 tvec(end)], 'FontSize', FontSize)
ylabel('FR (Hz)', 'FontSize', FontSize)

% tnew = AHV_tsd.tvec - AHV_tsd.tvec(1);
% AHV_tsdnew = tsd(tnew, AHV_tsd.data);
yyaxis right
plot(AHV_tsd.tvec, AHV_tsd.data)
ylabel('AHV (deg./sec)')

events_ts = LoadEvents([]);
% endtime = events_ts.t{2}(end) - events_ts.t{1}(1);
c = axis;
% axis([c(1) endtime c(3) c(4)]);

%% #8 MultiRaster
plot8 = subplot(4,6,13:18);

cfg = [];
cfg.uint = '64';
spikefiles = FindFiles('*.t');
cfg.fc = {spikefiles{iCell}};
Sp = LoadSpikes(cfg);

title(strcat(b, c), 'FontSize', FontSize)
cfg_mr = [];
cfg_mr.lfp = AHV_tsd;
cfg_mr.openNewFig = 0;
h = MultiRaster(cfg_mr, S);  % this is a hack
% h = MultiRaster(cfg_mr, Sp);
set(gca, 'FontSize', FontSize)
ylabel('AHV')
set(gca, 'YTickLabel', [])
c = axis;
% axis([c(1) endtime c(3) c(4)]);

%% #9 AHV and horizontal eye position
plot9 = subplot(4,6,19:24);
SSN = HD_GetSSN;
if exist(strcat(SSN, '-VT1_proc.mat'))
    load(strcat(SSN, '-VT1_proc.mat'), 'pupil');         % load the output of facemap
    % pupiltime = [1:length(pupil{1}.area)] ./ VT1fps;   % this is slightly off. Need to load timestamps from the Nvt file
    
    events_ts = LoadEvents([]);
    assert(strcmp(events_ts.label{1}, 'Starting Recording'))=1;  % assert is not working. matlab thinks its a variable 
%     starttime = events_ts.t{1}(1);
    
%     [~, videofn, ext] = fileparts(FindFiles('*VT1.nvt'));
%     cfg = [];
%     cfg.fn = strcat(videofn, ext);
%     cfg.removeZeros = 0 ;
%     pos_tsd = LoadPos(cfg);
    %     pupiltime = pos_tsd.tvec;   % it apprears to Nvt file is 2 frames longer than the number of frames from facemap
    %     pupiltime = pupiltime - pupiltime(1);
    
    [a, b, c] = fileparts(FindFile('*VT1.smi'));
    fn = strcat(b,c);
    tvec_raw = read_smi(fn);
    tvec = tvec_raw - starttime;
    
    yyaxis left
    plot(tvec, pupil{1}.com(:,2));         % pupil{1}.com(:,2)  is the horizontal eye position from facemap
    ylabel('Horiz. Eye Pos. (pixels)', 'FontSize', FontSize)
    xlabel('Time (sec)', 'FontSize', FontSize)
    set(gca, 'FontSize', FontSize)
    
    yyaxis right
    plot(AHV_tsd.tvec, AHV_tsd.data)
    ylabel('AHV (deg./sec)')
    c = axis;
    axis([c(1) endtime-starttime c(3) c(4)]);
    %%
    linkaxes([plot7 plot8 plot9], 'x');
else
end