function getSlowPhaseTC3(cfg_in, sd, iCell)
% JJS. 2024-03-12.
% Calculate and plot the AHV tuning curve
% input:   sd - session data structure with spike trains S and tsd of angular head velocity
% output:  tc_out - a structure
%               .usr.tc - nCells x mBins
%               .usr.occ_hist  -  occupancy in each bin
%               .usr.spk_hist  -  number of spikes in each AHV bin
%               .usr.good_idx  - bins with at least the minOcc # of samples
%               .usr.no_occ_idx  - bins with NaN values
%               .usr.binEdges
%               .usr.binCenters
%               .usr.nBins  - number of bins

% HELP
% Set 'doPlot'=1 to plot each TC. The default setting will smooth the data.
% One parameter that you may want to change is minOCC. This value determines the minimum number of samples (each 5 ms) needed in a given AHV bin to include that data.
% minOcc of 100 requires only 0.5 s of occupancy to include those values. I usually set this to 200 to include at least 1 s of data. But one may want to set it higher.
% JJS. 2024-03-15. This version restricts the spike train first and then does all operations with that restricted version. 


cfg_def = [];
cfg_def.doPlot = 1;
cfg_def.smooth = 1;
cfg_def.nBins = 100;
cfg_def.binEdges = {linspace(-200, 200, 101)};
cfg_def.occ_dt = median(diff(sd.AHV.tvec));
cfg_def.minOcc = 200;  % remember that Occ is measured in samples, not in seconds. Usually 5ms per sample, b/c the platform encoder sampling rate is 200Hz,
cfg_def.saccade_pre = .2;
cfg_def.saccade_post = .1;
cfg_def.LineWidth = 3;
cfg_def.FontSize = 20;
cfg_def.insetText = 18;
cfg_def.tightX = .025;
cfg_def.tightY = .035;
cfg_out = ProcessConfig2(cfg_def, cfg_in);

if exist(strcat(sd.SSN, '-saccades-edited.mat')) == 2
    load(strcat(sd.SSN, '-saccades-edited.mat')) % this will load the edited (curated) saccade timestamps and the pupil trace time series (position and velocity)
else
    error('saccade mat file not found')
end
sd = LoadSessionData([]);
myCell = SelectTS([], sd.S, iCell); % myCell.t{1}(keep{1}==1)=NaN;   % not sure if this line works 

timestouse = ~isnan(combinedSaccades);
combinedSaccadesToUse = combinedSaccades(timestouse);
[in, keep] = restrict(myCell, combinedSaccadesToUse - cfg_out.saccade_pre, combinedSaccadesToUse + cfg_out.saccade_post);  
keepers = keep{1,1} == 0;  % invert the output so that logical values of 1 indicate data that was not restricted (i.e. everything but peri-saccade times)
myCellr = myCell;
myCellr.t{1,1}(keep{1,1}==1) = NaN;

%% calculate AHV Tuning Curve
cfg_tc = [];
cfg_tc.nBins = 100;
cfg_tc.binEdges = {linspace(-200, 200, 101)};
cfg_tc.occ_dt = median(diff(AHV_tsd.tvec));
cfg_tc.minOcc = 100;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
tc_outAHV = TuningCurves(cfg_tc, myCellr, AHV_tsd);

%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #1 plot AHV tuning curve
% calculate raw firing rates 
clf; hold on;
p = subtightplot(4,2,1, [cfg_out.tightX cfg_out.tightY]);
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = cfg_tc.occ_dt;
cfg_Q.tvec_edges = AHV_tsd.tvec(1) : cfg_Q.dt : AHV_tsd.tvec(end);
F = MakeQfromS(cfg_Q, myCellr); % convert to FR
F.data = F.data ./ cfg_Q.dt;
F_idx = nearest_idx3(AHV_tsd.tvec, F.tvec);
AHV_F = F.data(:,F_idx);
ymax = max(AHV_F);
set(gca, 'TickDir', 'out', 'FontSize', cfg_out.FontSize)
plot(AHV_tsd.data, AHV_F, '.', 'MarkerSize', .5); hold on
set(gca, 'Ylim', [0 ymax], 'FontSize', cfg_out.FontSize)

% Add Tuning Curve
if cfg_out.smooth
    plot(tc_outAHV.binCenters, smoothdata(tc_outAHV.tc), 'LineWidth', cfg_out.LineWidth, 'Color', 'k');
else
    plot(tc_outAHV.binCenters, tc_outAHV.tc, 'LineWidth', 3, 'Color', 'k');
end
xlabel('AHV (deg/s)', 'FontSize', cfg_out.FontSize)
ylabel('FR (Hz)', 'FontSize', cfg_out.FontSize)
set(groot, 'DefaultLegendInterpreter', 'none')
title('Slow Phase Only')
text(NaN, NaN, 'CW', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.25 .85 0])
text(NaN, NaN, 'CCW', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.6 .85 0])
p.XAxisLocation = 'top';
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')

%% #2 Eye Velocity Tuning Curve
% Calculate pupil velocity firing rates
p = subtightplot(4,2,2, [cfg_out.tightX cfg_out.tightY]);
hold on
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = sd.tsdH_dt;
cfg_Q.tvec_edges = diffH.tvec(1): cfg_Q.dt: diffH.tvec(end); %#ok<*NODEF>
F = MakeQfromS(cfg_Q, myCellr); % convert to FR
F.data = F.data ./ cfg_Q.dt;
% find FR corresponding to each pupil position sample
F_idx = nearest_idx3(diffH.tvec, F.tvec);
tsdH_F = F.data(:,F_idx);
plot(diffH.data, tsdH_F(1,:), '.', 'MarkerSize', 2);   %  'color', [.8 .8 .8]
if size(diffH.tvec) == size(diffH.data)
    diffH.tvec = diffH.tvec'; % change the shape so that it is a "well-formed tsd" for tuning curves
end
set(gca, 'FontSize', cfg_out.FontSize)
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
text(-30, c(4)/2, 'nasal', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.15 .85 0])
text(30, c(4)/2, 'temporal', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.75 .85 0])
p.XAxisLocation = 'top';
xlabel('Eye Velocity (pixels/s)', 'FontSize', cfg_out.FontSize)
ylabel('FR (Hz)', 'FontSize', cfg_out.FontSize)
set(gca, 'FontSize', cfg_out.FontSize)
title(sd.SSN)
clear F

% calculate pupil velocity Tuning Curve
cfg_V = [];
cfg_V.nBins = 100;
cfg_V.binEdges = {linspace(-35, 35, 51)};
cfg_V.occ_dt = median(diff(diffH.tvec));
cfg_V.minOcc = 50;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
tc_vel = TuningCurves(cfg_V, myCellr, diffH);
% plot the tuning curve
if cfg_out.smooth
    plot(smoothdata(tc_vel.binCenters), tc_vel.tc, 'LineWidth', 3, 'Color', 'k');
else
    plot(tc_vel.binCenters, tc_vel.tc, 'LineWidth', 3, 'Color', 'k');
end
R = max(tc_vel.tc);
c = axis;
axis([-25 25 c(3) c(4)])
% c = axis;
% axis([-20 20 c(3) R + 10]); c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')

%% #3 Head Direction Tuning Curve
% Calculate head direction firing rates
% subplot(2,2,3)
[~, orientation, samplingrate, dt] = GetOrientationValues([]);
p = subtightplot(4,2,3, [cfg_out.tightX cfg_out.tightY]);
% calculate raw firing rates 
hold on;
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = cfg_tc.occ_dt;
cfg_Q.tvec_edges = orientation.tvec(1) : cfg_Q.dt : orientation.tvec(end);
F = MakeQfromS(cfg_Q, myCellr); % convert to FR
F.data = F.data ./ cfg_Q.dt;
F_idx = nearest_idx3(orientation.tvec, F.tvec);
AHV_F = F.data(:,F_idx);
ymax = max(AHV_F);
set(gca, 'TickDir', 'out', 'FontSize', cfg_out.FontSize)
plot(orientation.data, AHV_F, '.', 'MarkerSize', .5); hold on
set(gca, 'Ylim', [0 ymax], 'FontSize', cfg_out.FontSize)
% calculate head direction Tuning Curve
cfg_tc = [];
cfg_tc.nBins = 50;
cfg_tc.binEdges = {linspace(-200, 200, 51)};
cfg_tc.occ_dt = median(diff(orientation.tvec));
cfg_tc.minOcc = 100;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
tc_outHD = TuningCurves(cfg_tc, myCellr, orientation);
if cfg_out.smooth
    plot(tc_outHD.binCenters, smoothdata(tc_outHD.tc), 'LineWidth', cfg_out.LineWidth, 'Color', 'k');
else
    plot(tc_outHD.binCenters, tc_outHD.tc, 'LineWidth', 3, 'Color', 'k');
end
xlabel('Head Direction (deg)', 'FontSize', cfg_out.FontSize)
ylabel('FR (Hz)', 'FontSize', cfg_out.FontSize)
set(groot, 'DefaultLegendInterpreter', 'none')
% title('Slow Phase Only')
set(gca, 'FontSize', cfg_out.FontSize)
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')

%% #4 Eye Position Tuning Curve
% Calculate eye position firing rates
% subplot(2,2,4)
p = subtightplot(4,2,4, [cfg_out.tightX cfg_out.tightY]);
% calculate raw firing rates 
hold on;
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = cfg_tc.occ_dt;
cfg_Q.tvec_edges = tsdH.tvec(1) : cfg_Q.dt : tsdH.tvec(end);
F = MakeQfromS(cfg_Q, myCellr); % convert to FR
F.data = F.data ./ cfg_Q.dt;
F_idx = nearest_idx3(tsdH.tvec, F.tvec);
AHV_F = F.data(:,F_idx);
ymax = max(AHV_F);
set(gca, 'TickDir', 'out', 'FontSize', cfg_out.FontSize)
plot(tsdH.data, AHV_F, '.', 'MarkerSize', .5); hold on
set(gca, 'Ylim', [0 ymax], 'FontSize', cfg_out.FontSize)
% Add Tuning Curve
cfg_tc = [];
cfg_tc.nBins = 20;
cfg_tc.binEdges = {linspace(-40, 40, 21)};
cfg_tc.occ_dt = median(diff(tsdH.tvec));
cfg_tc.minOcc = 50;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
tsdH.data = tsdH.data'; % tuning curves wants a 1 x n, n x 1 size tsd
tc_outAHV = TuningCurves(cfg_tc, myCellr, tsdH);
if cfg_out.smooth
    plot(tc_outAHV.binCenters, smoothdata(tc_outAHV.tc), 'LineWidth', cfg_out.LineWidth, 'Color', 'k');
else
    plot(tc_outAHV.binCenters, tc_outAHV.tc, 'LineWidth', 3, 'Color', 'k');
end
xlabel('Eye Position (pixels)', 'FontSize', cfg_out.FontSize)
ylabel('FR (Hz)', 'FontSize', cfg_out.FontSize)
set(groot, 'DefaultLegendInterpreter', 'none')
text(NaN, NaN, 'nasal', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.25 .85 0])
text(NaN, NaN, 'temporal', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.6 .85 0])
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
axis tight
clear keepers
p = subtightplot(4,2,[5 7], [cfg_out.tightX cfg_out.tightY]);
%% #5 SpikePETH 
[outputS, ~, ~, outputIT, cfg_spikes] = SpikePETHvdm([], myCellr, combinedSaccadesToUse);
set(gca, 'FontSize', cfg_out.FontSize)

p = subtightplot(4,2,[6 8], [cfg_out.tightX cfg_out.tightY]);
plotSaccadeTrace([]);
title(sd.fn{iCell})

end




