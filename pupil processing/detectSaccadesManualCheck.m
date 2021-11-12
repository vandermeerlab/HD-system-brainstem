function [temporalSaccades, nasalSaccades, combinedSaccades, tsdH, tsdV, diffH, diffV, XT, YT, XN, YN, cfg] = detectSaccadesManualCheck(cfg_in)
% 2021-11. JJS.
% This function calculates saccades times from the eye position trace, similar to processPupilData2.m. Here, the threshold to use is manually adjusted,
%      and the trace is scrolled through to check/add/remove indivudal saccades that automatic method may have missed.

% INPUTS: 
%           cfg_in: define variables such as the number of pixels 
% OUTPUTS:
%           temporalSaccades:   timestamps for temporal saccades
%           nasalSaccades:      timestamps for nasal saccades
%           combinedSaccades:   both timestamps combined, sorted
%           tsdH:               tsd of horizontal pupil position
%           tsdV:               tsd of vertical pupil position
%           diffH:              tsd of horizontal pupil velocity  (*figure out units here) 
%           diffV:              tsd of vertical pupil velocity 
%           XT, YT:             x and y values for manually added TEMPORAL saccades 
%           XN, YN:             x and y values for manually added NASAL saccades 
%           cfg:                record of parameters used for analysis, like thresholds 

SSN = HD_GetSSN;
FontSize = 20;
cfg_def = [];
cfg_def.threshAdj  = 4;  % how many timesteps around a saccade that are disqualified from being considered as subsequent saccades. Saccades usually have a rebound that can hit the other threshold.
cfg_def.threshT = 10;  % positive displacement in image pixel space. TEMPORAL saccades.
cfg_def.threshN = -10; % negative displacement in image pixel space. NASAL saccades.

cfg_def.scalingfactor = 1;  % for shrinking the pupil trace so its the same height as diffH
cfg_def.artifactThresh = 4;  % units of pixels
cfg_def.doPlotThresholds = 1;
cfg_def.doPlotEverything = 1;
cfg = ProcessConfig2(cfg_def,cfg_in);  % not sure if there is a function difference btwn ver2 and ver1

%% Get timestamps for the pupil trace
% [~, videofn, ext] = fileparts(FindFiles('*VT1.nvt'));
% cfg_video = [];
% cfg_video.fn = strcat(videofn, ext);
% cfg_video.removeZeros = 0 ;
%
% pos_tsd = LoadPos(cfg_video);
% pupiltime = pos_tsd.tvec;   % it apprears that the Nvt file is 2 frames longer than the number of frames from facemap
%% Load Events and get the session start time
events_ts = LoadEvents([]);
index = strfind(events_ts.label, 'Starting Recording');
if index{1} == 1                                 % Start Recording should be in the first or second .label position.
    starttime = events_ts.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
elseif index{2} == 1
    starttime = events_ts.t{2}(1); % for session with laser events, the start recording time moves to position 2.
else
    error('could not find start time for this session')
end
%% Get timestamps from the .smi file
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
diffH.cfg.hdr{1}.Fs = 1 / median(diff(diffH.tvec));   % append the sampling rate

tstart = diffH.tvec(1);
tend = diffH.tvec(end);

%% Remove Artifacts from Camera Movement (vertical displacement)
% There shouldn't be much displacement in the vertical direction. Visual inspection showed that times with high power btwn 10-14 Hz were indicative of camera jitter that were not in fact saccades.
cfg.low_freq = 10;
cfg.high_freq = 14;
eeg1= diffV.data;
samp_freq = diffH.cfg.hdr{1}.Fs;
order = round(samp_freq); %determines the order of the filter used
if mod(order,2)~= 0
    order = order-1;
end
Nyquist=floor(samp_freq/2);%determines nyquist frequency
MyFilt=fir1(order,[cfg.low_freq cfg.high_freq]/Nyquist); %creates filter
filtered1 = Filter0(MyFilt,eeg1); %filters eeg1 between low_freq and high_freq
filt_hilb1 = hilbert(filtered1); %calculates the Hilbert transform of eeg1
amp1 = abs(filt_hilb1);%calculates the instantaneous amplitude of eeg1 filtered between low_freq and high_freq
amp1=amp1-mean(amp1); %removes mean of the signal because the DC component of a signal does not change the correlation
artifactIndex = amp1 > cfg.artifactThresh;  % find timepoints where power is high (suspect times for movement artifact)
suspectPoints = find(artifactIndex);

%% Thresholding the TEMPORAL saccades (positive AHV segments)
tP = diffH.data > cfg.threshT;  % data points above threshold
[~, index_tP] = find(tP);  % tvec indices for data points above threshold
[val_discard_tP, ~] = intersect(index_tP, suspectPoints);
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
temporalAmplitudes = diffH.datat(index_tP_final); 
nasalSaccades = diffH.tvec(index_nP_final);
nasalAmplitudes = diffH.data(index_nP_final); 
combinedSaccades = sort(horzcat(temporalSaccades, nasalSaccades));

disp(strcat('Num opposite peaks = ', num2str(sum(oppPeaks))));
disp(strcat('Adjacent opposite points removed  = ', num2str(length(indexes))));
disp(strcat('Num temporal saccades = ', num2str(length(index_tP_final))));
disp(strcat('Num nasal saccades = ', num2str(length(index_nP_final))));
disp(strcat('Total num saccades = ', num2str(length(index_nP_final)+length(index_tP_final))));

%% Plot the data and manually inspect
clf;
hold on
plot(diffH.tvec, diffH.data)
plot(diffV.tvec, diffV.data, 'm')
hold on
xlabel('Time (sec)', 'FontSize', FontSize)
ylabel('diff pupil pos', 'FontSize', FontSize)
title(SSN)
line([tstart tend], [cfg.threshT cfg.threshT], 'Color', 'k')
line([tstart tend], [cfg.threshN cfg.threshN], 'Color', 'k')
plot(diffH.tvec(index_tP_final), diffH.data(index_tP_final), 'r.', 'MarkerSize', 25)
plot(diffH.tvec(index_nP_final), diffH.data(index_nP_final), 'g.', 'MarkerSize', 25)
set(gca, 'FontSize', FontSize)

disp('find extra TEMPORAL saccades') 
count = 0;
XT= [];
YT = [];
while(1)
    fprintf(1, '\n');
    m = input('Do you want to continue, y/n?','s');
    if m == 'n'
        break
    end
    count = count + 1;
    fprintf(1, '\n');
    disp('Zoom into region of interest. Press return to continue')
    zoom on; 
    pause() % you can zoom with your mouse and when your image is okay, you press any key
    zoom off; % to escape the zoom mode
    fprintf(1, '\n');
    disp('Select point(s) for missing TEMPORAL saccade. Press return when finished.') 
    [x,y] =ginput;
    plot(x, y, 'k.', 'MarkerSize', 25)
    XT(end+1:end+length(x)) = x; 
    YT(end+1:end+length(y)) = y; 
    disp('point selected')
    clear m 
end

disp('find extra NASAL saccades') 
count = 0;
XN= [];
YN = [];
while(1)
    fprintf(1, '\n');
    m = input('Do you want to continue, y/n?','s');
    if m == 'n'
        break
    end
    count = count + 1;
    fprintf(1, '\n');
    disp('Zoom into region of interest. Press return to continue')
    zoom on; 
    pause() % you can zoom with your mouse and when your image is okay, you press any key
    zoom off; % to escape the zoom mode
    fprintf(1, '\n');
    disp('Select point(s) for missing TEMPORAL saccade. Press return when finished.') 
    [x,y] =ginput;
    plot(x, y, 'k.', 'MarkerSize', 25)
    XN(end+1:end+length(x)) = x; 
    YN(end+1:end+length(y)) = y; 
    disp('point selected')
    clear m 
end

temporalSaccades = sort(cat(1, temporalSaccades, XT)); 
nasalSaccades = sort(cat(1, nasalSaccades, XN)); 
combinedSaccades = sort(cat(1, temporalSaccades, nasalSaccades)); 




