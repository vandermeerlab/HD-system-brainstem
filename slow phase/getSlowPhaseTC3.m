function [numSpikesRemoved, data_out] = getSlowPhaseTC3(cfg_in, sd, iCell)
% JJS. 2024-03-12.
% Calculate and plot tuning curves for important variables during slow phase eye movements only. 
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
% 2024-03-22. Changed so that firing rate F is only computed once. Function now reports occ_dt, how many samples are removed, etc. Checks mean firing rate for Q and TCs.
% 2024-03-27. Changed restrict so that I am restricting out the time between saccades and using the tsd output instead of restricting the saccade intervals and using the keep output. 

format short
cfg_def = [];
cfg_def.doPlot = 1;
cfg_def.smooth = 0;
cfg_def.nBins = 100;
cfg_def.binEdges = {linspace(-200, 200, 101)};
cfg_def.occ_dt = median(diff(sd.AHV.tvec));
cfg_def.minOcc = 200;  % remember that Occ is measured in samples, not in seconds. Usually 5ms per sample, b/c the platform encoder sampling rate is 200Hz,
cfg_def.saccade_pre = .2;  % how many seconds to cut out before the saccade.
cfg_def.saccade_post = .1; % how many seconds to cut out after the saccade. 
cfg_def.LineWidth = 3;
cfg_def.FontSize = 20;
cfg_def.insetText = 18;
cfg_def.tightX = .075;
cfg_def.tightY = .045;

cfg_Q = [];
cfg_Q.dt = 0.05; % binsize in s
cfg_Q.smooth = 'gauss'; % [], 'gauss'
cfg_Q.gausswin_size = 1; % gaussian window size in seconds; only used if cfg.smooth = 'gauss'
cfg_Q.gausswin_sd = 0.1; % SD for gaussian convolution in seconds
cfg_def.cfg_Q = cfg_Q;

cfg_out = ProcessConfig2(cfg_def, cfg_in);

if exist(strcat(sd.SSN, '-saccades-edited.mat'),'file') == 2
    load(strcat(sd.SSN, '-saccades-edited.mat')) % this will load the edited (curated) saccade timestamps and the pupil trace time series (position and velocity)
else
    error('saccade mat file not found')
end
sd = LoadSessionData([]); % initialize the variables for this session
myCell = SelectTS([], sd.S, iCell); % choose a single neuron to work with
meanFRoverall = length(myCell.t{1,1})/sd.SessLength; 

%% Limit the data to everything but the quick phase periods 
timestouse = ~isnan(combinedSaccades); % create a logical with not NaN saccade times
combinedSaccadesToUse = combinedSaccades(timestouse); % select only non-NaN saccade values to use

pre = combinedSaccadesToUse - cfg_out.saccade_pre; 
post = combinedSaccadesToUse + cfg_out.saccade_post; 

[~, keep] = restrict(myCell, pre, post);  % get the indices for times around the moment of saccade

myCellr = myCell;
myCellr.t{1,1}(keep{1,1}==1) = [];  % remove the values that are within the saccade window 

disp(strcat('Total SPIKES = ' , num2str(length(sd.S.t{iCell}))))
numSpikesRemoved = sum((keep{1,1}==1)); disp(strcat('num SPIKES removed =', num2str(numSpikesRemoved)))
disp(strcat('Fraction = ', num2str(numSpikesRemoved/length(sd.S.t{iCell}),3)))
fprintf(1, '\n');

%% also output restricted firing rate (useful for GLM later)
Q = MakeQfromS(cfg_out.cfg_Q, myCellr);

[~, keep] = restrict(Q, pre, post);  
keep = (keep == 0);  % invert the output

fr = Q; fr.tvec = fr.tvec(keep); fr.data = fr.data(keep) ./ cfg_out.cfg_Q.dt;
data_out.fr = fr;

%% AHV
tStartAHV = horzcat(0, post);
tEndAHV = horzcat(pre, sd.AHV.tvec(end)); 
sizeFullAHV = length(sd.AHV.tvec);
[sacc,~] = restrict(sd.AHV, pre, post);  % this is the tsd of restricted saccade times 
size_sacc_AHV = length(sacc.tvec);
[ISI_AHV,~] = restrict(sd.AHV, tStartAHV, tEndAHV);  % ISI = inter-saccade interval 
size_ISI_AHV = length(ISI_AHV.tvec); 
total = size_sacc_AHV + size_ISI_AHV; 
assert(total == sizeFullAHV)
AHV_samplingrate = 1/median(diff(ISI_AHV.tvec)); disp(strcat('AHV sampling rate = ', num2str(AHV_samplingrate)))
fprintf(1, '\n');
AHVsamplesRemovedr = sum(sizeFullAHV - length(ISI_AHV.tvec)); 

[~,keepAHV] = restrict(sd.AHV, pre, post);  % AHV
keepersAHV = keepAHV == 0;  % invert the output so that logical values of 1 indicate data that was not restricted (i.e. everything but peri-saccade times)
AHVr = sd.AHV;
AHVr.tvec = sd.AHV.tvec(keepersAHV);
AHVr.data = sd.AHV.data(keepersAHV);
AHV_samplingrate = 1/median(diff(AHVr.tvec)); disp(strcat('AHV sampling rate = ', num2str(AHV_samplingrate)))
fprintf(1, '\n');
AHVsamplesRemoved = sum(keepAHV); 

%% EYE VELOCITY
[~,keepEV] = restrict(diffH, pre, post);  
keepersEV = keepEV == 0;  % invert the output so that logical values of 1 indicate data that was not restricted (i.e. everything but peri-saccade times)
diffHr = diffH;
diffHr.tvec = diffH.tvec(keepersEV);
diffHr.data = diffH.data(keepersEV);
diffH_samplingrate = 1/median(diff(diffHr.tvec)); disp(strcat('Eye velocity sampling rate = ', num2str(diffH_samplingrate)))
fprintf(1, '\n');
EVsamplesRemoved = sum(keepEV); 
%% HEAD DIRECTION
[~,keepHD] = restrict(sd.orientation, pre, post);  
keepersHD = keepHD == 0;  % invert the output so that logical values of 1 indicate data that was not restricted (i.e. everything but peri-saccade times)
orientationR = sd.orientation;
orientationR.tvec = orientationR.tvec(keepersHD);
orientationR.data = orientationR.data(keepersHD);
disp(strcat('Head direction sampling rate = ', num2str(sd.orientationsamplingrate))) 
fprintf(1, '\n');
HDsamplesRemoved = sum(keepHD); 
%% EYE POSITION
if size(tsdH.tvec) == size(tsdH.data) %#ok<NODEF>
    tsdH.data = tsdH.data';
end
[~,keepEP] = restrict(tsdH, pre, post);   
keepersEP = keepEP == 0;  % invert the output so that logical values of 1 indicate data that was not restricted (i.e. everything but peri-saccade times)
tsdHr = tsdH;
tsdHr.tvec = tsdHr.tvec(keepersEP);
tsdHr.data = tsdHr.data(keepersEP);
tsdH_samplingrate = 1/median(diff(tsdH.tvec)); disp(strcat('Eye position sampling rate = ', num2str(tsdH_samplingrate)))
fprintf(1, '\n');
EPsamplesRemoved = sum(keepEP); 

disp(strcat('AHV samples Removed = ', num2str(AHVsamplesRemoved)))
disp(strcat('HD samples Removed = ', num2str(HDsamplesRemoved)))
fprintf(1, '\n');
disp(strcat('EV samples Removed = ', num2str(EVsamplesRemoved)))
disp(strcat('EP samples Removed = ', num2str(EPsamplesRemoved)))
fprintf(1, '\n');

data_out.horiz_eye_pos = tsdHr;
data_out.horiz_eye_vel = diffHr;
data_out.AHV = AHVr;
data_out.S = myCellr;

%% get AHV Tuning Curve
cfg_tcAHV = [];
cfg_tcAHV.nBins = 100;
cfg_tcAHV.binEdges = {linspace(-200, 200, 101)};
cfg_tcAHV.minOcc = 100;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
cfg_tcAHV.occ_dt = median(diff(AHVr.tvec));
tc_outAHV = TuningCurves(cfg_tcAHV, myCellr, AHVr);

%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #1 AHV tuning curve
% calculate raw firing rates 
clf; hold on;
% subplot(2,2,1)
subtightplot(4,2,1, [cfg_out.tightX cfg_out.tightY]);
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = cfg_tcAHV.occ_dt;
cfg_Q.tvec_edges = AHVr.tvec(1) : cfg_Q.dt : AHVr.tvec(end);
F = MakeQfromS(cfg_Q, myCell); % convert to FR
F.data = F.data ./ cfg_Q.dt;
F_idx = nearest_idx3(AHVr.tvec, F.tvec);
AHV_F = F.data(:,F_idx);
ymax = max(AHV_F);
set(gca, 'TickDir', 'out', 'FontSize', cfg_out.FontSize)
plot(AHVr.data, AHV_F, '.', 'MarkerSize', .5); hold on
set(gca, 'Ylim', [0 ymax], 'FontSize', cfg_out.FontSize)

% add the Tuning Curve
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
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')

%% #2 Eye Velocity Tuning Curve
% Calculate pupil velocity firing rates
subtightplot(4,2,2, [cfg_out.tightX cfg_out.tightY]);
hold on
% find FR corresponding to each pupil position sample
F_idxEV = nearest_idx3(diffHr.tvec, F.tvec);
EV_F = F.data(:,F_idxEV);
ymaxEV = max(EV_F);
plot(diffHr.data, EV_F(1,:), '.', 'MarkerSize', 2);   %  'color', [.8 .8 .8]
if size(diffHr.tvec) == size(diffHr.data)
    diffHr.tvec = diffHr.tvec'; % change the shape so that it is a "well-formed tsd" for tuning curves
end

set(gca, 'FontSize', cfg_out.FontSize)
% title('Pupil Position (pixels)')
% axis tight
c = axis;
axis([-45 45 c(3) c(4)])
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
text(-30, c(4)/2, 'nasal', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.15 .85 0])
text(30, c(4)/2, 'temporal', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.75 .85 0])
xlabel('Eye Velocity (pixels/s)', 'FontSize', cfg_out.FontSize)
ylabel('FR (Hz)', 'FontSize', cfg_out.FontSize)
set(gca, 'FontSize', cfg_out.FontSize)
title(sd.fn{iCell,1})

% add the Tuning Curve
cfg_V = [];
cfg_V.nBins = 100;
cfg_V.binEdges = {linspace(-25, 25, 100)};
cfg_V.occ_dt = median(diff(diffHr.tvec)); 
cfg_V.minOcc = 25;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
tc_velEV = TuningCurves(cfg_V, myCellr, diffHr);

if cfg_out.smooth
    plot(smoothdata(tc_velEV.binCenters), tc_velEV.tc, 'LineWidth', 3, 'Color', 'k');
else
    plot(tc_velEV.binCenters, tc_velEV.tc, 'LineWidth', 3, 'Color', 'k');
end
c = axis;
axis([-10 10 c(3) max(F.data)]); c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')

%% #3 Head Direction Tuning Curve
% Calculate head direction firing rates
% subplot(2,2,3)
subtightplot(4,2,3, [cfg_out.tightX cfg_out.tightY]);
% calculate raw firing rates 
hold on;
F_idxHD = nearest_idx3(orientationR.tvec, F.tvec);
HD_F = F.data(:,F_idxHD);
ymaxHD = max(HD_F);
set(gca, 'TickDir', 'out', 'FontSize', cfg_out.FontSize)
plot(orientationR.data, HD_F, '.', 'MarkerSize', .5); hold on
set(gca, 'Ylim', [0 ymaxHD], 'FontSize', cfg_out.FontSize)

% add the Tuning Curve
cfg_tcHD = [];
cfg_tcHD.nBins = 50;
cfg_tcHD.binEdges = {linspace(-200, 200, 51)};
cfg_tcHD.occ_dt = median(diff(orientationR.tvec)); 
cfg_tcHD.minOcc = 100;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
tc_outHD = TuningCurves(cfg_tcHD, myCellr, orientationR);
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
subtightplot(4,2,4, [cfg_out.tightX cfg_out.tightY]);
% calculate raw firing rates 
hold on;
F_idxEP = nearest_idx3(tsdHr.tvec, F.tvec);
EP_F = F.data(:,F_idxEP);
ymaxEP = max(EP_F);
set(gca, 'TickDir', 'out', 'FontSize', cfg_out.FontSize)
plot(tsdHr.data, EP_F, '.', 'MarkerSize', .5); hold on
set(gca, 'Ylim', [0 ymaxEP], 'FontSize', cfg_out.FontSize)

% Add Tuning Curve
cfg_tcEP = [];
cfg_tcEP.nBins = 20;
cfg_tcEP.binEdges = {linspace(-40, 40, 21)};
cfg_tcEP.occ_dt = median(diff(tsdHr.tvec)); 
cfg_tcEP.minOcc = 50;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
tc_outEP = TuningCurves(cfg_tcEP, myCellr, tsdHr);
if cfg_out.smooth
    plot(tc_outEP.binCenters, smoothdata(tc_outEP.tc), 'LineWidth', cfg_out.LineWidth, 'Color', 'k');
else
    plot(tc_outEP.binCenters, tc_outEP.tc, 'LineWidth', 3, 'Color', 'k');
end
xlabel('Eye Position (pixels)', 'FontSize', cfg_out.FontSize)
ylabel('FR (Hz)', 'FontSize', cfg_out.FontSize)
set(groot, 'DefaultLegendInterpreter', 'none')
text(NaN, NaN, 'nasal', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.25 .85 0])
text(NaN, NaN, 'temporal', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.6 .85 0])
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
axis tight

%% #5  Plot the spikePETH 
subtightplot(4,2,[5 7], [cfg_out.tightX cfg_out.tightY]);
cfg_in.doPlot = 1;
cfg_in.doBar = 0;
[~, ~, ~, ~, ~] = SpikePETHvdm([], myCellr, combinedSaccadesToUse);
set(gca, 'FontSize', cfg_out.FontSize)
xlabel('time peri Saccade (s)')


%% Display occ_dt
disp(strcat('occ_dt for AHV = ', num2str(cfg_tcAHV.occ_dt,3)))
disp(strcat('occ_dt for HD = ', num2str(cfg_tcHD.occ_dt,3)))
disp(strcat('occ_dt for EV = ', num2str(cfg_V.occ_dt,3)))
disp(strcat('occ_dt for EP = ', num2str(cfg_tcEP.occ_dt,3)))
fprintf(1, '\n');

%% #6 Plot the Firing Rate Scatterplot and Tuning Curve for Unrestricted AHV data 
% calculate raw firing rates 
subtightplot(4,2,[6 8], [cfg_out.tightX cfg_out.tightY]); hold on
cfg_Qu = [];
cfg_Qu.smooth = 'gauss';
cfg_Qu.gausswin_sd = 0.05;
cfg_Qu.dt = median(diff(AHV_tsd.tvec));
cfg_Qu.tvec_edges = AHV_tsd.tvec(1) : cfg_Q.dt : AHV_tsd.tvec(end);
Fu = MakeQfromS(cfg_Qu, myCell); % convert to FR
Fu.data = Fu.data ./ cfg_Qu.dt;
F_idxu = nearest_idx3(AHV_tsd.tvec, Fu.tvec);
AHV_Fu = F.data(:,F_idxu);
ymaxu = max(AHV_Fu);
set(gca, 'TickDir', 'out', 'FontSize', cfg_out.FontSize)
plot(sd.AHV.data, AHV_Fu, '.', 'MarkerSize', .5); hold on
set(gca, 'Ylim', [0 ymaxu], 'FontSize', cfg_out.FontSize)

% get AHV Tuning Curve
cfg_tcAHVu = [];
cfg_tcAHVu.nBins = 100;
cfg_tcAHVu.binEdges = {linspace(-200, 200, 101)};
cfg_tcAHVu.minOcc = 100;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
cfg_tcAHVu.occ_dt = median(diff(sd.AHV.tvec));
tc_outAHVu = TuningCurves(cfg_tcAHVu, myCell, sd.AHV);
% plot the Tuning Curve
if cfg_out.smooth
    plot(tc_outAHVu.binCenters, smoothdata(tc_outAHVu.tc), 'LineWidth', cfg_out.LineWidth, 'Color', 'k');
else
    plot(tc_outAHVu.binCenters, tc_outAHVu.tc, 'LineWidth', 3, 'Color', 'k');
end
xlabel('AHV (deg/s)', 'FontSize', cfg_out.FontSize)
ylabel('FR (Hz)', 'FontSize', cfg_out.FontSize)
set(groot, 'DefaultLegendInterpreter', 'none')
text(NaN, NaN, 'CW', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.25 .85 0])
text(NaN, NaN, 'CCW', 'FontSize', cfg_out.insetText, 'Units', 'normalized', 'Position', [.6 .85 0])
c = axis;
line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
text(150, 20, 'All Data', 'FontSize', 20)

%% Display the mean firing rates
fprintf(1, '\n');
disp(strcat('Overall FR = ', num2str(meanFRoverall,4), '_Hz')) % note: change this to sprintf
disp(strcat('Q matrix mean F = ', num2str(nanmean(F.data),4), '_Hz')) % note: change this to sprintf

fprintf(1, '\n');
meanAHV_F = nanmean(AHV_F);
meanAHV_TC = nansum(tc_outAHV.occ_hist .* tc_outAHV.tc) ./ nansum(tc_outAHV.occ_hist); % average firing rate, weighted by occupancy
disp(strcat('mean AHV_F = ', num2str(meanAHV_F,4), '_Hz')) % note: change this to sprintf
disp(strcat('mean AHV tc = ', num2str(meanAHV_TC,4), '_Hz')) % note: change this to sprintf
subtightplot(4,2,1, [cfg_out.tightX cfg_out.tightY]);
c = axis;
line([c(1) c(2)], [meanAHV_F meanAHV_F], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
line([c(1) c(2)], [meanAHV_TC meanAHV_TC], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
fprintf(1, '\n');

meanEV_F = mean(EV_F);
meanEV_TC = nansum(tc_velEV.occ_hist .* tc_velEV.tc) ./ nansum(tc_velEV.occ_hist);  % average firing rate, weighted by occupancy
disp(strcat('mean EV_F = ', num2str(meanEV_F,4), '_Hz')) % note: change this to sprintf
disp(strcat('mean EV tc = ', num2str(meanEV_TC,4), '_Hz')) % note: change this to sprintf
subtightplot(4,2,2, [cfg_out.tightX cfg_out.tightY]);
c = axis;
line([c(1) c(2)], [meanEV_F meanEV_F], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
line([c(1) c(2)], [meanEV_TC meanEV_TC], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
fprintf(1, '\n');

meanHD_F = nanmean(HD_F);
meanHD_TC = nansum(tc_outHD.occ_hist .* tc_outHD.tc) ./ nansum(tc_outHD.occ_hist);
disp(strcat('mean HD_F = ', num2str(meanHD_F,4), '_Hz')) % note: change this to sprintf
disp(strcat('mean HD tc = ', num2str(meanHD_TC,4), '_Hz')) % note: change this to sprintf
subtightplot(4,2,3, [cfg_out.tightX cfg_out.tightY]);
c = axis;
line([c(1) c(2)], [meanHD_F meanHD_F], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
line([c(1) c(2)], [meanHD_TC meanHD_TC], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
fprintf(1, '\n');

meanEP_F = nanmean(EP_F);
meanEP_TC = nansum(tc_outEP.occ_hist .* tc_outEP.tc) ./ nansum(tc_outEP.occ_hist);
disp(strcat('mean EP_F = ', num2str(meanEP_F,4), '_Hz')) % note: change this to sprintf
disp(strcat('mean EP tc = ', num2str(meanEP_TC,4), '_Hz')) % note: change this to sprintf
subtightplot(4,2,4, [cfg_out.tightX cfg_out.tightY]);
c = axis;
line([c(1) c(2)], [meanEP_F meanEP_F], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
line([c(1) c(2)], [meanEP_TC meanEP_TC], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')
fprintf(1, '\n');

meanAHV_Fu = nanmean(AHV_Fu);
meanAHVu_TC = nansum(tc_outAHVu.occ_hist .* tc_outAHVu.tc) ./ nansum(tc_outAHVu.occ_hist);

disp(strcat('mean AHV all FR = ', num2str(meanAHV_Fu,4), '_Hz')) % note: change this to sprintf
disp(strcat('mean AHV all tc = ', num2str(meanAHVu_TC,4), '_Hz')) % note: change this to sprintf
subtightplot(4,2,[6 8], [cfg_out.tightX cfg_out.tightY]);
c = axis;
line([c(1) c(2)], [meanAHV_Fu meanAHV_Fu], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--')
line([c(1) c(2)], [meanAHVu_TC meanAHVu_TC], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--')


end

