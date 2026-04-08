function detectSaccadesManualCheck3(cfg_in)
% 2021-11. JJS.
% This function calculates saccades times from the eye position trace, similar to processPupilData2.feedback. Here, the threshold to use is manually adjusted,
%      and the trace is scrolled through to check/add/remove indivudal saccades that automatic method may have missed.
% This version accepts all potential saccades in the first pass, and sorts them (temporal vs. nasal) afterward, according to sign (+ = temp., - = nasal).
% *** Need to also add a way to remove saccades that are incorrect
%% THIS  VERSION CANNOT APPEND NEW CHANGES TO AN EXISTING FILE!!!!

% INPUTS:
%           cfg_in: define variables such as the number of pixels
% OUTPUTS:
%           m.temporalSaccades:   timestamps for temporal saccades
%           m.nasalSaccades:      timestamps for nasal saccades
%           m.combinedSaccades:   both timestamps combined, sorted
%           m.tsdH:               tsd of horizontal pupil position
%           m.tsdV:               tsd of vertical pupil position
%           m.diffH:              tsd of horizontal pupil velocity  (*figure out units here)
%           m.diffV:              tsd of vertical pupil velocity
%           m.XS, YS:             x and y values for manually added TEMPORAL saccades
%           m.XN, YN:             x and y values for manually added NASAL saccades
%           m.temporalAmplitdues: amplitude of temporal saccades
%           m.nasalAmplitudes:    amplitude of nasal saccades
%           m.cfg:                record of parameters used for analysis, like thresholds

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
cfg_def.LineWidth = 3;
cfg = ProcessConfig2(cfg_def,cfg_in);  % not sure if there is a function difference btwn ver2 and ver1

%% Get timestamps for the pupil trace
% [~, videofn, ext] = fileparts(FindFiles('*VT1.nvt'));
% cfg_video = [];
% cfg_video.fn = strcat(videofn, ext);
% cfg_video.removeZeros = 0 ;
%
% pos_tsd = LoadPos(cfg_video);
% pupiltime = pos_tsd.tvec;   % it apprears that the Nvt file is 2 frames longer than the number of frames from facemap
%% Calculate AHV
cfg_AHV = [];
[AHV_tsd] = Get_AHV(cfg_AHV);

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
cfg.high_freq = 15;
eeg1= diffV.data;
samp_freq = diffH.cfg.hdr{1}.Fs;
order = round(samp_freq); %determines the order of the filter used
if mod(order,2)~= 0
    order = order-1;
end
Nyquist=floor(samp_freq/2);%determines nyquist frequency
if Nyquist == 15
    disp('Frame Rate is 30 Hz')
    cfg.high_freq = 14; % this is a hacky. A few sessions were run with a frame rate of 30Hz. Matlab won't let me filter at or above 15 Hz for these sessions.
end
MyFilt=fir1(order,[cfg.low_freq cfg.high_freq]/Nyquist); %creates filter
filtered1 = Filter0(MyFilt,eeg1); %filters eeg1 between low_freq and high_freq
filt_hilb1 = hilbert(filtered1); %calculates the Hilbert transform of eeg1
amp1 = abs(filt_hilb1);%calculates the instantaneous amplitude of eeg1 filtered between low_freq and high_freq
amp1tsd  = tsd(diffV.tvec, amp1);
amp1=amp1-mean(amp1); %removes mean of the signal because the DC component of a signal does not change the correlation
artifactIndex = amp1 > cfg.artifactThresh;  % find timepoints where power is high (suspect times for movement artifact)
suspectPoints = find(artifactIndex);

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
MyFilt=fir1(order,[cfg.low_freq cfg.high_freq]/Nyquist); %creates filter
filtered2 = Filter0(MyFilt,eeg2); %filters eeg1 between low_freq and high_freq
%% This is an addition on 2021/11/12. hilbert.m gives all NaN values if thery are ANY NaNs in the input.
F2 = fillmissing(filtered2,'constant', 0);
filt_hilb2 = hilbert(F2); %calculates the Hilbert transform of eeg1. Note *** Hilbert transform cannot accept NaNs. Use interp to fill these in.
amp2 = abs(filt_hilb2);%calculates the instantaneous amplitude of eeg1 filtered between low_freq and high_freq
amp2=amp2-mean(amp2); %removes mean of the signal because the DC component of a signal does not change the correlation
amp2tsd  = tsd(diffH.tvec, amp2);

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
temporalAmplitudes = diffH.data(index_tP_final);
nasalSaccades = diffH.tvec(index_nP_final);
nasalAmplitudes = diffH.data(index_nP_final);

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
% legend('horiz eye vel.', 'vertical eye vel.', 'horizontal eye position', 'filtered vert. vel. 10-15 Hz', 'filtered horiz. vel. 10-15 Hz', '', '', '')
yyaxis right
plot(AHV_tsd.tvec, AHV_tsd.data, 'Color', [.75 .75 0])
% c = axis;
% line([c(1) c(2)], [0 0], 'Color', [.75 .75 0], 'LineStyle', '--', 'LineWidth', 3)        % plotting another line makes it glitchy with the
ylabel('horizontal pupil position', 'FontSize', FontSize)
yyaxis left

fprintf(1, '\n');
fprintf(1, '\n');
disp('find extra saccades of any type')
count = 0;
XS= [];
YS = [];
Button = [];
while(1)
    fprintf(1, '\n');
    feedback = input('Do you want to continue, y/n? ...','s');
    if feedback == 'n'
        break
    end
    clear feedback
    count = count + 1;
    fprintf(1, '\n');
    disp('Zoom into region of interest. Press return to stop zoom.')
    zoom on;
    pause() % you can zoom with your mouse and when your image is okay, you press any key
    zoom off; % to escape the zoom mode
    fprintf(1, '\n');
    disp('Select point(s) for missing saccade(s). Use [r] key to indicate REMOVAL. Press return when finished.')
    [x,y, button] =ginput;  % button will be a number, according to ASCII format. Lower case 'r' is the number [114].
    removalIndex = button == 114;
    plot(x, y, 'k.', 'MarkerSize', 25)
    text(x(removalIndex), y(removalIndex), 'R', 'FontSize', 20)
    accept = input('Accept these selections?  enter/n ...', 's');
    if strcmp(accept,'') || isempty(accept)
        XS(end+1:end+length(x)) = x;
        YS(end+1:end+length(y)) = y;
        if ~isempty(Button)
            Button(end+1:end+length(button)) = button;
        else
            Button = button;
        end
        disp('point selected')
    else
        disp('Continue through prompts to re-select points')
        plot(x, y, 'w.', 'MarkerSize', 25)
    end
end
%% Sort out the button presses to see if there were any saccade REMOVAL entries
RemovalIndex = Button == 114;   % Need to remove saccades for this index that are closest to the chose XY,YS values.
XS_Remove = XS(RemovalIndex);
YS_Remove = YS(RemovalIndex);

Removeflag = 1;
if ~isempty(YS_Remove)
    Ytemporal_Remove = YS_Remove > 0;
    Ynasal_Remove = YS_Remove < 0;
    
    assert(sum(Ytemporal_Remove) + sum(Ynasal_Remove) == length(YS_Remove))
else
    Removeflag = 0;          % no saccades to REMOVE
end
%% Remove values
iCount = 0;
if Removeflag == 1
    if sum(Ytemporal_Remove) > 0
        for iC = 1:length(Ytemporal_Remove)
            iCount = iCount + 1;
            [minValueT,closestIndexT] = min(abs(XS_Remove(iC)-temporalSaccades));
            XS_Remove_updated(iCount) = temporalSaccades(closestIndexT);
            YS_Remove_updated(iCount) = temporalAmplitudes(closestIndexT);
            temporalSaccades(closestIndexT) = NaN;
            temporalAmplitudes(closestIndexT) = NaN;
        end
    end
    if sum(Ynasal_Remove) > 0
        for iC = 1:length(Ynasal_Remove)
            iCount = iCount + 1;
            [minValueN,closestIndexN] = min(abs(XS_Remove(iC)-nasalSaccades));
            XS_Remove_updated(iCount) = nasalSaccades(closestIndexN);
            YS_Remove_updated(iCount) = nasalAmplitudes(closestIndexN);
            nasalSaccades(closestIndexN) = NaN;
            nasalAmplitudes(closestIndexN) = NaN;
        end
    end
end
%% Add Values
XS_Add = XS(~RemovalIndex);  % x coordinates for points that were ADDED
YS_Add = YS(~RemovalIndex);

Addflag = 1;
if ~isempty(YS_Add)
    Ytemporal_Add = YS_Add > 0;      % index of which YS (or XS) values to add
    Ynasal_Add = YS_Add < 0;
    assert(sum(Ytemporal_Add) + sum(Ynasal_Add) == length(YS_Add))
else
    Addflag = 0;             % no saccades to ADD
end
if Addflag == 1
    [temporalSaccades_sorted, sortT] = sort(cat(2, temporalSaccades, XS_Add(Ytemporal_Add)));
    [nasalSaccades_sorted, sortN] = sort(cat(2, nasalSaccades, XS_Add(Ynasal_Add)));
    
    temporalAmplitudes_temp = cat(2, temporalAmplitudes, YS_Add(Ytemporal_Add));
    nasalAmplitudes_temp = cat(2, nasalAmplitudes, YS_Add(Ynasal_Add));
    
    temporalAmplitudes_sorted =  temporalAmplitudes_temp(sortT);
    nasalAmplitudes_sorted =  nasalAmplitudes_temp(sortN);
else
    temporalSaccades_sorted = temporalSaccades; temporalAmplitudes_sorted = temporalAmplitudes;
    nasalSaccades_sorted = nasalSaccades; nasalAmplitudes_sorted = nasalAmplitudes;
end
%% Plot the data and manually inspect
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
plot(temporalSaccades_sorted, temporalAmplitudes_sorted, 'r.', 'MarkerSize', 25)
plot(nasalSaccades_sorted, nasalAmplitudes_sorted, 'g.', 'MarkerSize', 25)
% Plot the points that were added
if Addflag == 1
    plot(XS_Add(Ytemporal_Add), YS_Add(Ytemporal_Add), 'ro', 'MarkerSize', 35, 'LineWidth', 5)   % orange cross
    plot(XS_Add(Ynasal_Add), YS_Add(Ynasal_Add), 'go', 'MarkerSize', 35, 'LineWidth', 5)    % forest green cross
end
% plot the points that were removed
if Removeflag == 1
    plot(XS_Remove_updated, YS_Remove_updated, 'ko', 'MarkerSize', 35, 'LineWidth', 5)   % black open circles
end
set(gca, 'FontSize', FontSize)
% legend('horiz eye vel.', 'vertical eye vel.', 'horizontal eye position', 'filtered vert. vel. 10-15 Hz', 'filtered horiz. vel. 10-15 Hz', '', '', '')
yyaxis right
plot(AHV_tsd.tvec, AHV_tsd.data, 'Color', [.75 .75 0])
ylabel('horizontal pupil position', 'FontSize', FontSize)
yyaxis left

temporalSaccades = temporalSaccades_sorted;
nasalSaccades = nasalSaccades_sorted;
temporalAmplitudes = temporalAmplitudes_sorted;
nasalAmplitudes = nasalAmplitudes_sorted;
num_temporalSaccades = length(~isnan(temporalSaccades));
num_nasalSaccades = length(~isnan(nasalSaccades_sorted));

disp('Inspect the output. Press any key to continue') 
pause
%% Put data into structure
% m.temporalSaccades = temporalSaccades_sorted;
% m.nasalSaccades = nasalSaccades_sorted;
% m.combinedSaccades = combinedSaccades;
% m.temporalAmplitudes = temporalAmplitudes_sorted;
% m.nasalAmplitudes = nasalAmplitudes_sorted;
% m.tsdH = tsdH;
% m.tsdV = tsdV;
% m.diffH = diffH;
% m.diffV = diffV;
% m.XS_Remove = XS_Remove;   % x coordinates for saccades that were REMOVED
% m.YS_Remove = YS_Remove;
% m.XS_Add = XS_Add;         % x coordinates for saccades that were ADDED
% m.YS_Add = YS_Add;
% m.cfg = cfg;
% m.tstart = tstart;
% m.tend = tend;
% m.tvec = tvec;
% m.Button = Button;     % ASCII codes for button presses / mouse clicks to select saccades
%% Save data: overwrite

disp('Saving data as new -saccade.m file')
save(strcat(SSN, '-saccades.mat'), 'temporalSaccades', 'num_temporalSaccades', 'nasalSaccades', 'num_nasalSaccades', 'temporalAmplitudes', 'nasalAmplitudes', 'tsdH', 'tsdV', 'diffH', 'diffV', 'XS_Remove', 'YS_Remove', 'XS_Add', 'YS_Add', 'cfg', 'tstart', 'tend', 'tvec', 'Button', 'amp1tsd', 'amp2tsd', 'index_tP_final', 'index_nP_final', 'AHV_tsd');
savefig(strcat(SSN, '-saccades.fig'))


