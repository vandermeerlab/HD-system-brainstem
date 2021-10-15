function plot_AHVtuningMegaplot(iCell, varargin)

hist2Xmax = 0.1;
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
S = LoadSpikes(cfg);
if subtractStartTime == 1 % New cheetah versions have timestamps 
    events_ts = LoadEvents([]);
    assert(strcmp(events_ts.label{1}, 'Starting Recording'))=1;
    for iC = 1:length(S.t) 
        S.t{iC} = S.t{iC} - events_ts.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time. 
    end
end
% get AHV
cfg_AHV = [];
cfg_AHV.subsample_factor = 10;
[AHV_tsd, tc_out] = AHV_tuning(S, cfg_AHV);
AHV_dt = median(diff(AHV_tsd.tvec));

%% #2 plot scatterplot
subplot(3,6,2)
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = AHV_dt;
F = MakeQfromS(cfg_Q, S); % convert to FR
% convert to FR
F.data = F.data ./ cfg_Q.dt;
% find FR corresponding to each AHV sample
F_idx = nearest_idx3(AHV_tsd.tvec, F.tvec);
AHV_F = F.data(:,F_idx);
ymax = max(AHV_F(iCell,:));

set(gca, 'TickDir', 'out', 'FontSize', FontSize)
plot(AHV_tsd.data, AHV_F(iCell,:), '.', 'MarkerSize', .5);
set(gca, 'Ylim', [0 ymax], 'FontSize', FontSize)
xlabel('AHV (deg./sec)', 'FontSize', FontSize)
ylabel('FR (Hz)', 'FontSize', FontSize)
title('Scatterplot')
h = get(gca, 'XLim');
text(.75*h(1), 10, 'CW', 'FontSize', 12)
text(.5*h(2), 10, 'CCW', 'FontSize', 12)

%% #1 plot tuning curves
subplot(3,6,1)
fc = FindFiles('*.t', 'CheckSubdirs', 0);
[a, b, c] = fileparts(fc{iCell});
plot(tc_out.usr.binCenters, tc_out.tc(iCell,:), 'LineWidth', LineWidth);
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
subplot(3,6,3)
cfg_acf = [];
cfg_acf.binsize = 0.001;
cfg_acf.max_t = 0.5;
cfg_acf.smooth = 1; % set to 1 to compute ccf on SDF, 0 for raw spikes
cfg_acf.gauss_w = 1; % width of Gaussian convolution window (in s)
cfg_acf.gauss_sd = 0.005; % SD of Gaussian convolution window (in s)
set(gca, 'TickDir', 'out', 'FontSize', FontSize)
[acf, tvec] = ccf(cfg_acf, S.t{iCell}, S.t{iCell});
midpoint = ceil(length(acf)./2);
acf(midpoint) = 0;
plot(tvec, acf, 'LineWidth', 1);
xlabel('Time (sec)', 'FontSize', FontSize)
set(gca, 'xtick', [-.5 -.25 0 .25 .5], 'FontSize', FontSize); grid on;
title('Acorr')


%% #4 acf zoomed in 
subplot(3,6,4); hold on
cfg_acf = [];
cfg_acf.binsize = 0.001;
cfg_acf.max_t = 0.05;
cfg_acf.smooth = 1; % set to 1 to compute ccf on SDF, 0 for raw spikes
cfg_acf.gauss_w = 1; % width of Gaussian convolution window (in s)
cfg_acf.gauss_sd = 0.005; % SD of Gaussian convolution window (in s)
set(gca, 'TickDir', 'out', 'FontSize', FontSize)
[acf, tvec] = ccf(cfg_acf, S.t{iCell}, S.t{iCell});
midpoint = ceil(length(acf)./2);
acf(midpoint) = 0;
plot(tvec, acf, 'LineWidth', 1);
xlabel('Time (sec)', 'FontSize', FontSize)
set(gca, 'xtick', [-.05 0 .05], 'FontSize', FontSize); grid on;
title('Acorr')


%% #5 HistISI
subplot(3,6,5); hold on
[h, n] = HistISIsubplot(S.t{iCell});
HistISIsubplot(S.t{iCell});
c = axis;
[~, i] = max(h);
line([n(i) n(i)], [0 c(4)], 'color', 'k');
grid on
set(gca, 'TickDir', 'out', 'XLim', [histXmin histXmax], 'FontSize', FontSize)
xlabel('Time (sec)', 'FontSize', FontSize)
title('HistISI')
% ylabel('Count')


%% #6 HistISI zoomed in 
subplot(3,6,6); hold on






%%  #7 Firing Rate
plot7 = subplot(3,6,7:12); hold on
cfg_Q = []; cfg_Q.dt = 0.001; cfg_Q.gausswin_sd = 0.05;cfg_Q.smooth = 'gauss';
Q = MakeQfromS(cfg_Q, S);
tvec = Q.tvec - Q.tvec(1);
yyaxis left 
plot(tvec, Q.data(iCell,:)./cfg_Q.dt)
maxFR = max(Q.data(iCell,:)./cfg_Q.dt);
set(gca, 'Xlim', [0 tvec(end)], 'FontSize', FontSize)
ylabel('FR (Hz)', 'FontSize', FontSize)

% tnew = AHV_tsd.tvec - AHV_tsd.tvec(1);
% AHV_tsdnew = tsd(tnew, AHV_tsd.data);
yyaxis right
plot(AHV_tsd.tvec, AHV_tsd.data)
ylabel('AHV (deg./sec)')


%% #8 MultiRaster
plot8 = subplot(3,6,13:18)
title(strcat(b, c), 'FontSize', FontSize)
cfg_mr = [];
cfg_mr.lfp = AHV_tsd;
cfg_mr.openNewFig = 0;
h = MultiRaster(cfg_mr, S);
set(gca, 'FontSize', FontSize)
ylabel('AHV')
set(gca, 'YTickLabel', [])

linkaxes([plot7 plot8], 'x')
