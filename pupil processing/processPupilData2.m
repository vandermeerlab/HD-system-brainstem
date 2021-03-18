function [temporalSaccades, nasalSaccades, tsdH, tsdV, diffH, diffV, filteredH, filteredV] = processPupilData2(cfg_in)

% JJS. 2021-03-13.
% Basic working version that finds candidate saccade events and filters out camera jitter detected by vertical displacement

% This version takes the diff first and then filters
SSN = HD_GetSSN;
FontSize = 20;
cfg_def = [];
cfg_def.threshH = 10;
cfg_def.threshL = -10;
cfg_def.doPlotThresholds = 1;
cfg_def.doPlotEverything = 0;
cfg_def.removeArtifcat = 1;
cfg_def.scalingfactor = 1;  % for shrinking the pupil trace so its the same height as diffH
cfg_def.artifactThresh = 4;  % units of pixels sq.
cfg = ProcessConfig(cfg_def,cfg_in);

%% Get AHV trace for later plotting
[~, orientation, ~, ~] = GetOrientationValues([]);
[orientationtouse] = downsampleOrientationValues(orientation, 10);
[AHVtsd] = GetAHV_values(orientationtouse);

%% Get timestamps for the pupil trace
[~, videofn, ext] = fileparts(FindFiles('*VT1.nvt'));
cfg_video = [];
cfg_video.fn = strcat(videofn, ext);
cfg_video.removeZeros = 0 ;

pos_tsd = LoadPos(cfg_video);
pupiltime = pos_tsd.tvec;   % it apprears that the Nvt file is 2 frames longer than the number of frames from facemap
pupiltime = pupiltime - pupiltime(1);
pupiltime = pupiltime(2:end-1); % this is maybe a temporary fix. Need to find out why there is a 2 frames discrepancy. See how many frames in mp4

%% Get the pupil trace dervied from Facemap
f = FindFiles('*VT1_proc.mat');
load(f{1}, 'pupil');
% Subtract the mean
meanH = mean(pupil{1}.com(:,2));
meanV = mean(pupil{1}.com(:,1));
% Make it into a TSD
tsdH = tsd(pupiltime, pupil{1}.com(:,2) - meanH);   % tsd of horizontal pupil position
tsdV = tsd(pupiltime, pupil{1}.com(:,1) - meanV);   % tsd of vertical pupil position

diffH = tsd(pupiltime(2:end), diff(tsdH.data)');     % Should this be (1:end-1) or (2:end)?
diffV = tsd(pupiltime(2:end), diff(tsdV.data)');

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
if cfg.removeArtifcat == 1
cfg.low_freq = 10;
cfg.high_freq = 15;
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
diff_tP = horzcat([NaN diff(index_tP)]);  % spacing for thresholded data
adjacent_tP = diff_tP ==1;       % which values of index_tP are adjacent (one timestep apart)
index_tP_temp = index_tP(adjacent_tP == 0); 
temporalSaccades = diffH.tvec(index_tP_temp);
[val_discard_tP, ~] = intersect(index_tP_temp, suspectPoints);
index_tP_touse = setdiff(index_tP_temp, val_discard_tP);
disp(strcat('Num nasal saccades = ', num2str(length(index_tP_touse))));

%% Thresholding the NASAL saccades (negative AHV segments)    
nP = diffH.data < cfg.threshL;  % data points above threshold
[~, index_nP] = find(nP);  % tvec indices for data points above threshold
diff_nP = horzcat([NaN diff(index_nP)]);  % spacing for thresholded data
adjacent_nP = diff_nP ==1;       % which values of index_tP are adjacent (one timestep apart)
index_nP_temp = index_nP(adjacent_nP == 0); 
nasalSaccades = diffH.tvec(index_nP_temp);
[val_discard_nP, ~] = intersect(index_nP_temp, suspectPoints);
index_nP_touse = setdiff(index_nP_temp, val_discard_nP);
disp(strcat('Num nasal saccades = ', num2str(length(index_nP_touse))));








%% Plot the data in a subplot
if cfg_def.doPlotEverything == 1
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
    linkaxes([ax ay az], 'xy')
end

if cfg_def.doPlotThresholds == 1
    clf
    hold on
    plot(diffH.tvec, diffH.data)
    plot(diffV.tvec, diffV.data, 'm')
    hold on
    xlabel('Time (sec)', 'FontSize', FontSize)
    ylabel('diff pupil pos', 'FontSize', FontSize)
    title(SSN)
    line([tstart tend], [cfg.threshH cfg.threshH], 'Color', 'k')
    line([tstart tend], [cfg.threshL cfg.threshL], 'Color', 'k')
    plot(diffH.tvec(index_tP_touse), diffH.data(index_tP_touse), 'r.', 'MarkerSize', 25)
    plot(diffH.tvec(index_nP_touse), diffH.data(index_nP_touse), 'g.', 'MarkerSize', 25)
    set(gca, 'FontSize', FontSize)
    
end


end



% suspectPoints = find(artifactIndex);
% tP = diffH.data > cfg.threshH;   % data points above threshold
% [~, index_tP] = find(tP);  % tvec indices for data points above threshold
% [val_discard_tP, ~] = intersect(index_tP, suspectPoints);  % identify overlap with suspect points
% index_tP_temp = setdiff(index_tP, val_discard_tP);   % remove indices of suspect points 
% diff_tP = horzcat([NaN diff(index_tP_temp)]);  % spacing for thresholded data
% adjacent_tP = diff_tP ==1;       % which values of index_tP are adjacent (one timestep apart)
% %% For adjacent points, choose the peak value
% 
% for iAdj = ind_adj 
%     temp = diffH.data(index_tP_temp([iAdj iAdj-1]));
%     [m, i] = max(temp);
%     if i == 2
%         index_tP_temp(iAdj-1) = NaN;
%     elseif i ==1
%         index_tP_temp(iAdj) = NaN; 
%     else
%         warning('peak for adjacent points not found')
%     end
% end


