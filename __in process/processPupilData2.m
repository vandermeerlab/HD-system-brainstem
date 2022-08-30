function [temporalSaccades, nasalSaccades, combinedSaccades, index_tP_final, index_nP_final, tsdH, tsdV, diffH, diffV, temporalAmplitudes, nasalAmplitudes] = processPupilData2(cfg_in, varargin)
% JJS. 2021-03-13.
% Remove jitter first then thresholds.
% This is the in progress version.

% INPUTS:
%           cfg_in: define variables such as the number of pixels
% OUTPUTS:
%           temporalSaccades:   timestamps for temporal saccades
%           nasalSaccades:      timestamps for nasal saccades
%           combinedSaccades:   both timestamps combined, sorted
%           index_tP_final:     indices for temporal saccades, in terms of the pupil position CSC
%           index_nP_final:     indices for nasal saccades, in terms of the pupil position CSC
%           tsdH:               tsd of horizontal pupil position
%           tsdV:               tsd of vertical pupil position
%           diffH:              tsd of horizontal pupil velocity  (*figure out units here)
%           diffV:              tsd of vertical pupil velocity
doSave = 1;
SSN = HD_GetSSN;
FontSize = 20;
cfg_def = [];
cfg_def.LineWidth = 3;
cfg_def.threshAdj  = 4;  % how many timesteps around a saccade that are disqualified from being considered as subsequent saccades. Saccades usually have a rebound that can hit the other threshold.
cfg_def.threshT = 1;  % positive displacement in image pixel space. TEMPORAL saccades.
cfg_def.threshN = -1; % negative displacement in image pixel space. NASAL saccades.

cfg_def.scalingfactor = 1;  % for shrinking the pupil trace so its the same height as diffH
cfg_def.artifactThresh = 4;  % units of pixels
cfg_def.doPlotThresholds = 1;
cfg = ProcessConfig2(cfg_def,cfg_in);  % not sure if there is a function difference btwn ver2 and ver1

%% Get timestamps for the pupil trace
% [~, videofn, ext] = fileparts(FindFiles('*VT1.nvt'));
% cfg_video = [];
% cfg_video.fn = strcat(videofn, ext);
% cfg_video.removeZeros = 0 ;
%
% pos_tsd = LoadPos(cfg_video);
% pupiltime = pos_tsd.tvec;   % it apprears that the Nvt file is 2 frames longer than the number of frames from facemap
events_ts = LoadEvents([]);
% sd = LoadSessionData([], 'EYE', false);
% assert(strcmp(events_ts.label{1}, 'Starting Recording')==1);
index = strfind(events_ts.label, 'Starting Recording');
if index{1} == 1                                 % Start Recording should be in the first or second .label position.
    starttime = events_ts.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
elseif index{2} == 1
    starttime = events_ts.t{2}(1); % for session with laser events, the start recording time moves to position 2.
else
    error('could not find start time for this session')
end

% [~, videofn, ext] = fileparts(FindFiles('*VT1.smi'));
% filename = strcat(videofn, ext);
% A = readtable(filename, 'FileType', 'text');   % grab the data from the .smi file, which contains the .mp4 timestamps
% B = table2array(A(:,6));                       %
% C = regexp(B,'[\d*\.]*\d*','match');
% D = [];
% for iCell = 1:length(C)
%     D(iCell,1) = str2double(cell2mat(C{iCell}));
% end
% E = D(1:end - 2);  % the last two values from the table are empty (not data) and turn into NANs. Discard those entries.
% pupiltime_raw = D*10^-6;    % convert from microseconds to seconds
% pupiltime = pupiltime_raw - starttime;

[~, b, c] = fileparts(FindFile('*VT1.smi'));
fn = strcat(b,c);
tvec_raw = read_smi(fn);
tvec = tvec_raw - starttime;

%% Get the pupil trace dervied from Facemap
f = FindFiles('*VT1_proc.mat');
load(f{1}, 'pupil');
pupilH = pupil{1}.com(:,2);
pupilV = pupil{1}.com(:,1);

% Subtract the mean
meanH = nanmean(pupilH);
meanV = nanmean(pupilV);
% Make it into a TSD
tsdH = tsd(tvec, pupilH - meanH);   % tsd of horizontal pupil position
tsdV = tsd(tvec, pupilV - meanV);   % tsd of vertical pupil position

diffH = tsd(tvec(2:end)', diff(tsdH.data)');     % Should this be (1:end-1) or (2:end)?
diffV = tsd(tvec(2:end)', diff(tsdV.data)');     % tsd of vertical pupil velocity

tstart = diffH.tvec(1);
tend = diffH.tvec(end);

%% Remove Artifacts from Camera Movement
% VERTICAL eye velocity
% There shouldn't be much displacement in the vertical direction. Visual inspection showed that times with high power btwn 10-14 Hz were indicative of camera jitter that were not in fact saccades.
cfg.low_freq = 10;
cfg.high_freq = 15;
eeg1= diffV.data;
samp_freq =  1 / median(diff(diffV.tvec));
order = round(samp_freq); %determines the order of the filter used
if mod(order,2)~= 0
    order = order-1;
end
Nyquist=floor(samp_freq/2);%determines nyquist frequency
if Nyquist == 15
    disp('Frame Rate is 30 Hz')
    cfg.high_freq = 14; % this is a hacky. A few sessions were run with a frame rate of 30Hz. Matlab won't let me filter at or above 15 Hz for these sessions.
end
%% Interp NaNs if any exist (otherwise filtering won't work)
eeg1new = eeg1;
nanx = isnan(eeg1);
t = 1:numel(eeg1);
eeg1new(nanx) = interp1(t(~nanx), eeg1(~nanx), t(nanx));

MyFilt=fir1(order,[cfg.low_freq cfg.high_freq]/Nyquist); %creates filter
filtered1 = Filter0(MyFilt,eeg1new); %filters eeg1 between low_freq and high_freq
%% This is an addition on 2021/11/12. hilbert.m gives all NaN values if thery are ANY NaNs in the input.
F = fillmissing(filtered1,'constant', 0);
filt_hilb1 = hilbert(F); %calculates the Hilbert transform of eeg1. Note *** Hilbert transform cannot accept NaNs. Use interp to fill these in.
amp1 = abs(filt_hilb1);%calculates the instantaneous amplitude of eeg1 filtered between low_freq and high_freq
amp1=amp1-mean(amp1); %removes mean of the signal because the DC component of a signal does not change the correlation
amp1tsd  = tsd(diffV.tvec, amp1);
artifactIndex = amp1 > cfg.artifactThresh;  % find timepoints where power is high (suspect times for movement artifact)
suspectPoints = find(artifactIndex);
formatSpec = 'number of VERTICAL points with suspect filtered power is is %4.0d points';
fprintf(formatSpec, length(suspectPoints));

%% Remove Artifacts from Camera Movement
% HORIZONTAL eye velocity
cfg.low_freq = 10;
cfg.high_freq = 15;
eeg2= diffH.data;
samp_freq =  1 / median(diff(diffV.tvec));
order = round(samp_freq); %determines the order of the filter used
if mod(order,2)~= 0
    order = order-1;
end
Nyquist=floor(samp_freq/2);%determines nyquist frequency
if Nyquist == 15
    disp('Frame Rate is 30 Hz')
    cfg.high_freq = 14; % this is a hacky. A few sessions were run with a frame rate of 30Hz. Matlab won't let me filter at or above 15 Hz for these sessions.
end
%% Interp NaNs if any exist (otherwise filtering won't work)
eeg2new = eeg2;
nanx = isnan(eeg2);
t = 1:numel(eeg2);
eeg2new(nanx) = interp1(t(~nanx), eeg2(~nanx), t(nanx));
MyFilt=fir1(order,[cfg.low_freq cfg.high_freq]/Nyquist); %creates filter
filtered2 = Filter0(MyFilt,eeg2new); %filters eeg1 between low_freq and high_freq
%% This is an addition on 2021/11/12. hilbert.m gives all NaN values if thery are ANY NaNs in the input.
F2 = fillmissing(filtered2,'constant', 0);
filt_hilb2 = hilbert(F2); %calculates the Hilbert transform of eeg1. Note *** Hilbert transform cannot accept NaNs. Use interp to fill these in.
amp2 = abs(filt_hilb2);%calculates the instantaneous amplitude of eeg1 filtered between low_freq and high_freq
amp2=amp2-mean(amp2); %removes mean of the signal because the DC component of a signal does not change the correlation
amp2tsd  = tsd(diffH.tvec, amp2);
artifactIndex2 = amp2 > cfg.artifactThresh;  % find timepoints where power is high (suspect times for movement artifact)
suspectPoints2 = find(artifactIndex2);
formatSpec = 'number of HORIZONTAL points with suspect filtered power is is %4.0d points';
fprintf(formatSpec, length(suspectPoints2));

%% Thresholding the TEMPORAL saccades (positive AHV segments)
tP = diffH.data > cfg.threshT;  % data points above threshold
[~, index_tP] = find(tP);  % tvec indices for data points above threshold
[val_discard_tP, ~] = intersect(index_tP, suspectPoints);
formatSpec = 'number of potential saccades removed due to filtering is %4.0d points';
fprintf(formatSpec, length(val_discard_tP));

index_tP_temp = setdiff(index_tP, val_discard_tP);
sac_amps_tP = diffH.data(index_tP_temp);
diff_tP = horzcat([NaN diff(index_tP_temp)]);  % spacing for thresholded data
adjacent_tP = diff_tP ==1;       % which values of index_tP are adjacent (one timestep apart)
new_adj_tP = adjacent_tP;
indices = find(adjacent_tP);
for iAdj = indices(1:end)
    [~, c] = max(sac_amps_tP(iAdj-1: iAdj));
    if c ==1
        new_adj_tP(iAdj-1) = 0;
        new_adj_tP(iAdj) = 1;
    elseif c == 2
        new_adj_tP(iAdj-1) = 1;
        new_adj_tP(iAdj) = 0;
    else
        fprintf(1, '\n');
        warning('cant find the max for saccade start')
    end
end
index_tP_touse = index_tP_temp(new_adj_tP == 0);

%% Thresholding the NASAL saccades (negative AHV segments)
nP = diffH.data < cfg.threshN;  % data points above threshold
[~, index_nP] = find(nP);  % tvec indices for data points above threshold
[val_discard_nP, ~] = intersect(index_nP, suspectPoints);
index_nP_temp = setdiff(index_nP, val_discard_nP);
sac_amps_nP = diffH.data(index_nP_temp);
diff_nP = horzcat([NaN diff(index_nP_temp)]);  % spacing for thresholded data
adjacent_nP = diff_nP ==1;       % which values of index_tP are adjacent (one timestep apart)
new_adj_nP = adjacent_nP;
indices = find(adjacent_nP);
for iAdj = indices(1:end)
    [~, c] = min(sac_amps_nP(iAdj-1: iAdj));
    if c ==1
        new_adj_nP(iAdj-1) = 0;    % ignore the first, keep the second
        new_adj_nP(iAdj) = 1;
    elseif c == 2
        new_adj_nP(iAdj-1) = 1;
        new_adj_nP(iAdj) = 0;
    else
        fprintf(1, '\n');
        warning('cant find the max for saccade start')
    end
end
index_nP_touse = index_nP_temp(new_adj_nP == 0);

%% Find Temproal and Nasal saccades that are within a few timesteps and remove
B = horzcat(index_tP_touse, index_nP_touse);   % combined array of temporal and nasal saccade indices
tPsize = 1:length(index_tP_touse);
nPsize = length(index_tP_touse)+1: length(index_tP_touse) + length(index_nP_touse);
[sortedB, I] = sort(B);                        % now in order of appearance
diffB = horzcat([NaN diff(sortedB)]);          % # of time steps between adjacent saccades
oppPeaks = diffB <= cfg.threshAdj;                         % find saccades that are 3 apart or less (60ms or less apart)
oppPeaksdata = diffH.data(sortedB);            % pupil positions for the sorted array
new_oppPeaks = oppPeaks;
idx = find(oppPeaks);
for iAdj = idx(1:end)
    [~, c] = min(abs(oppPeaksdata(iAdj-1: iAdj)));
    if c ==1
        new_oppPeaks(iAdj-1) = 1;
        new_oppPeaks(iAdj) = 0;
    elseif c == 2
        new_oppPeaks(iAdj-1) = 0;
        new_oppPeaks(iAdj) = 1;
    else
        fprintf(1, '\n');
        warning('cant find the max for saccade start')
    end
end
list = B(I(new_oppPeaks==1)) ;
c = ismember(B, list);
indexes = find(c);
B(indexes) = NaN;
tP_wnan = B(tPsize);
index_tP_final = tP_wnan(~isnan(tP_wnan));
nP_wnan = B(nPsize);
index_nP_final = nP_wnan(~isnan(nP_wnan));

temporalSaccades = diffH.tvec(index_tP_final);
temporalAmplitudes = diffH.data(index_tP_final);
nasalSaccades = diffH.tvec(index_nP_final);
nasalAmplitudes = diffH.data(index_nP_final);
combinedSaccades = sort(horzcat(temporalSaccades, nasalSaccades));

fprintf(1, '\n');
disp(strcat('Num opposite peaks = ', num2str(sum(oppPeaks))));
disp(strcat('Adjacent opposite points removed  = ', num2str(length(indexes))));
disp(strcat('Num temporal saccades = ', num2str(length(index_tP_final))));
disp(strcat('Num nasal saccades = ', num2str(length(index_nP_final))));
disp(strcat('Total num saccades = ', num2str(length(index_nP_final)+length(index_tP_final))));

if cfg.doPlotThresholds == 1
    clf;
    hold on
    plot(diffH.tvec, diffH.data)
    plot(diffV.tvec, diffV.data, 'm')
    plot(tsdH.tvec, tsdH.data, 'Color', [.301 .745 .933], 'LineStyle', '--')
    %     plot(tsdH.tvec, tsdH.data, 'Color', 'k', 'LineStyle', '--')
    plot(amp1tsd.tvec, amp1tsd.data, 'k', 'LineWidth', cfg.LineWidth)
    plot(amp2tsd.tvec, amp2tsd.data, 'Color', [.85 .325 .098], 'LineWidth', cfg.LineWidth)
    xlabel('Time (sec)', 'FontSize', FontSize)
    ylabel('diff pupil pos', 'FontSize', FontSize)
    title(SSN)
    line([tstart tend], [cfg_def.threshT cfg_def.threshT], 'Color', 'r')
    line([tstart tend], [cfg_def.threshN cfg_def.threshN], 'Color', 'g')
    line([tstart tend], [cfg.artifactThresh cfg.artifactThresh], 'Color', 'k', 'LineStyle', '--')
    line([tstart tend], [-cfg.artifactThresh -cfg.artifactThresh], 'Color', 'k', 'LineStyle', '--')
    plot(diffH.tvec(index_tP_final), diffH.data(index_tP_final), 'r.', 'MarkerSize', 25)
    plot(diffH.tvec(index_nP_final), diffH.data(index_nP_final), 'g.', 'MarkerSize', 25)
    set(gca, 'FontSize', FontSize)
    legend('horiz eye vel.', 'vertical eye vel.', 'horizontal eye position', 'filtered vert. vel. 10-15 Hz', 'filtered horiz. vel. 10-15 Hz')
end

if doSave == 1
    %% Save data
    disp('Saving data as new -saccade.mat file')
    save(strcat(SSN, '-saccades.mat'), 'temporalSaccades', 'nasalSaccades', 'combinedSaccades', 'index_tP_final', 'index_nP_final', 'tsdH', 'tsdV', 'diffH', 'diffH')
    if cfg.doPlotThresholds == 1
        savefig(strcat(SSN, '-saccades.fig'))
    end
end

%           temporalSaccades:   timestamps for temporal saccades
%           nasalSaccades:      timestamps for nasal saccades
%           combinedSaccades:   both timestamps combined, sorted
%           index_tP_final:     indices for temporal saccades, in terms of the pupil position CSC
%           index_nP_final:     indices for nasal saccades, in terms of the pupil position CSC
%           tsdH:               tsd of horizontal pupil position
%           tsdV:               tsd of vertical pupil position
%           diffH:              tsd of horizontal pupil velocity  (*figure out units here)
%           diffV:



% %% Filtering
% % The horizontal trace
% cfg_filter = [];
% cfg_filter.type = 'butter';
% cfg_filter.band = 'highpass';
% cfg_filter.f = 1;
% diffH.cfg.hdr{1}.Fs = 1 / median(diff(diffH.tvec));   % append the sampling rate
% filteredH = FilterLFP(cfg_filter, diffH);
%
% % The vertical trace
% cfg_filter = [];
% cfg_filter.type = 'butter';
% cfg_filter.band = 'highpass';
% cfg_filter.f = 1;
% diffV.cfg.hdr{1}.Fs = 1 / median(diff(diffV.tvec));   % append the sampling rate
% filteredV = FilterLFP(cfg_filter, diffV);

