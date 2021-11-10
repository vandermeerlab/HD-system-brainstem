function [temporalSaccades, nasalSaccades, combinedSaccades, index_tP_final, index_nP_final, tsdH, tsdV, diffH, diffV] = processPupilData2(cfg_in, sd, varargin)
% JJS. 2021-03-13.
% Remove jitter first then thresholds.
% This is the in progress version.

SSN = HD_GetSSN;
FontSize = 20;
cfg_def = [];
cfg_def.threshAdj  = 4;
cfg_def.threshH = 10;
cfg_def.threshL = -10;

cfg_def.scalingfactor = 1;  % for shrinking the pupil trace so its the same height as diffH
cfg_def.artifactThresh = 4;  % units of pixels sq.
cfg_def.doPlotThresholds = 1;
cfg_def.doPlotEverything = 1;
cfg = ProcessConfig(cfg_def,cfg_in);

%% Get timestamps for the pupil trace
% [~, videofn, ext] = fileparts(FindFiles('*VT1.nvt'));
% cfg_video = [];
% cfg_video.fn = strcat(videofn, ext);
% cfg_video.removeZeros = 0 ;
%
% pos_tsd = LoadPos(cfg_video);
% pupiltime = pos_tsd.tvec;   % it apprears that the Nvt file is 2 frames longer than the number of frames from facemap
% events_ts = LoadEvents([]);
% sd = LoadSessionData([], 'EYE', false);
% assert(strcmp(events_ts.label{1}, 'Starting Recording')==1);
index = strfind(sd.Events.label, 'Starting Recording');
if index{1} == 1                                 % Start Recording should be in the first or second .label position.
    starttime = sd.Events.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
elseif index{2} == 1
    starttime = sd.Events.t{2}(1);
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
diffV = tsd(tvec(2:end)', diff(tsdV.data)');

tstart = diffH.tvec(1);
tend = diffH.tvec(end);
%% Filtering
% The horizontal trace
cfg_filter = [];
cfg_filter.type = 'butter';
cfg_filter.band = 'highpass';
cfg_filter.f = 1;
diffH.cfg.hdr{1}.Fs = 1 / median(diff(diffH.tvec));   % append the sampling rate
filteredH = FilterLFP(cfg_filter, diffH);

% The vertical trace
cfg_filter = [];
cfg_filter.type = 'butter';
cfg_filter.band = 'highpass';
cfg_filter.f = 1;
diffV.cfg.hdr{1}.Fs = 1 / median(diff(diffV.tvec));   % append the sampling rate
filteredV = FilterLFP(cfg_filter, diffV);

%% Remove Artifacts from Camera Movement (vertical displacement)

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
artifactIndex = amp1 > cfg.artifactThresh;
suspectPoints = find(artifactIndex);

%% Thresholding the TEMPORAL saccades (positive AHV segments)
tP = diffH.data > cfg.threshH;  % data points above threshold
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
nP = diffH.data < cfg.threshL;  % data points above threshold
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
nasalSaccades = diffH.tvec(index_nP_final);
combinedSaccades = sort(horzcat(temporalSaccades, nasalSaccades));

disp(strcat('Num opposite peaks = ', num2str(sum(oppPeaks))));
disp(strcat('Adjacent opposite points removed  = ', num2str(length(indexes))));
disp(strcat('Num temporal saccades = ', num2str(length(index_tP_final))));
disp(strcat('Num nasal saccades = ', num2str(length(index_nP_final))));
disp(strcat('Total num saccades = ', num2str(length(index_nP_final)+length(index_tP_final))));
%% Plot the data in a subplot
if cfg.doPlotEverything == 1
    figure
    ax = subplot(3,1,1);
    plot(tsdH.tvec, tsdH.data./cfg.scalingfactor)
    ylabel('pupil pos.')
    title(SSN)
    
    ay = subplot(3,1,2);
    plot(diffH.tvec, diffH.data)
    ylabel('diff pupil pos')
    
    az = subplot(3,1,3);
    plot(filteredH.tvec, filteredH.data)
    ylabel('filtered diff pupil pos')
    xlabel('Time (sec)')
    
    %     linkaxes([ax ay], 'xy')
    linkaxes([ax ay], 'xy')
end

if cfg.doPlotThresholds == 1
    figure;
    hold on
    plot(diffH.tvec, diffH.data)
    plot(diffV.tvec, diffV.data, 'm')
    hold on
    xlabel('Time (sec)', 'FontSize', FontSize)
    ylabel('diff pupil pos', 'FontSize', FontSize)
    title(SSN)
    line([tstart tend], [cfg.threshH cfg.threshH], 'Color', 'k')
    line([tstart tend], [cfg.threshL cfg.threshL], 'Color', 'k')
    plot(diffH.tvec(index_tP_final), diffH.data(index_tP_final), 'r.', 'MarkerSize', 25)
    plot(diffH.tvec(index_nP_final), diffH.data(index_nP_final), 'g.', 'MarkerSize', 25)
    set(gca, 'FontSize', FontSize)
    
end



