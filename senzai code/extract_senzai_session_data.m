% 2024-10-24. JJS. Extracts Senzai data for one mouse.

cd('J:\senzai dataset\data files')
%% M77
%% Import the spike train data
disp('importing data')
M77_TCs = importdata('YutaTest77c_HeadDirection_OpenField.mat'); % data strcuture that contains Maps, kappa_ppln, pval_ppln, z_ppln.
% Maps contains the tuning curves. pval_ppln is presumably a p value for some significance test for head direction selectivity.
% Not sure what kappa is. z_ppln might be z-score.
numCellsM77 = length(M77_TCs.Maps);

spiketrain_1_M77 = importdata('YutaTest77c.res.1'); % Shank 1 from mouse 77. I think these spike trains are split in two because of their size.
spiketrain_1_M77 = spiketrain_1_M77./20000; % divide by the sampling rate (20kHz) to get seconds
spiketrain_2_M77 = importdata('YutaTest77c.res.2'); % Shank 2 from mouse 77.
spiketrain_2_M77 = spiketrain_2_M77./20000;

clu_1_M77 = importdata('YutaTest77c.clu.1');  % 'clu_1' here means the clusters from shank 1. 'clu_2' means the clusters from shank 2.
M77shank1_num = clu_1_M77(1); % this includes zeros and ones
clu_1_M77toUse = clu_1_M77(2:end);
clu_2_M77 = importdata('YutaTest77c.clu.2');
M77shank2_num = clu_2_M77(1); % this includes zeros and ones
clu_2_M77toUse = clu_2_M77(2:end);

assert(length(spiketrain_1_M77) == length(clu_1_M77toUse))
assert(length(spiketrain_2_M77) == length(clu_2_M77toUse))
%% Isolate the individual clusters (i.e., the single units)
clu_1_M77_IDs = unique(clu_1_M77toUse);
temp1 = ismember(clu_1_M77_IDs, 0:1); % 0 = artifact and 1 = noise. Remove these from consideration.
num_clu_1_M77 = sum(temp1==0); % how many single units (clusters) are there?
clu_1_M77_IDsToUse = clu_1_M77_IDs(clu_1_M77_IDs ~= 0 & clu_1_M77_IDs ~= 1); % Use clusters that are not zero or one.

clu_2_M77_IDs = unique(clu_2_M77toUse);
temp2 = ismember(clu_2_M77_IDs, 0:1); % 0 = artifact and 1 = noise. Remove these from consideration.
num_clu_2_M77 = sum(temp2==0); % how many single units (clusters) are there?
clu_2_M77_IDsToUse = clu_2_M77_IDs(clu_2_M77_IDs ~= 0 & clu_2_M77_IDs ~= 1); % Use clusters that are not zero or one.

totalCellsM77 = num_clu_1_M77 + num_clu_2_M77;
if numCellsM77 ~= totalCellsM77
    warning('M77 neurons do not add up')
end

%% Create the spike train structure
clear M77
counterToUse = 0;
M77.type = 'ts';
for iC = 1: totalCellsM77  % combine spike trains from each file into one structure
    if iC <= num_clu_1_M77
        M77.t{iC} = spiketrain_1_M77(clu_1_M77toUse == clu_1_M77_IDsToUse(iC));
    elseif iC > num_clu_1_M77 && iC <= totalCellsM77
        counterToUse = counterToUse + 1;
        M77.t{iC} = spiketrain_2_M77(clu_2_M77toUse == clu_2_M77_IDsToUse(counterToUse));
    else
        warning('problem with neuron count')
    end
    M77.label{iC} = 'M77';
end

%% Get heading and AHV
M77_heading = importdata('YutaTest77c_OpenField_HeadDirection.mat'); % contains Neck&NoseOmitIdx [value of 1 = omit], rho(distance between nose to neck),
% theta (heading, in radians), t (timestamp, IN SECONDS)
M77tracking_dur = (M77_heading.t(end) - M77_heading.t(1))/60; % How many minutes long the tracking was in the open field. This is 98 minutes for M77.
M77_Q = unwrap(M77_heading.theta);
M77_Q_deg = M77_Q*180/pi;
window = 0.1;
postsmoothing = .05;
M77_AHV = dxdt(M77_heading.t, M77_Q_deg, 'window', window, 'postsmoothing', postsmoothing);
M77_AHVtsd = tsd(M77_heading.t, M77_AHV);
M77_heading_sampling_rate = 1/median(diff(M77_heading.t));  % should be 50Hz
% Isolate low AHV times
AHVthresh = 10; % cm/sec
M77_low_AHV = abs(M77_AHV) < AHVthresh;
M77_low_AHV_samples = sum(M77_low_AHV);
M77_low_AHV_seconds = M77_low_AHV_samples/M77_heading_sampling_rate;
M77_low_AHV_minutes = M77_low_AHV_seconds/60;
fraction_M77_low_AHV = M77_low_AHV_samples / length(M77_heading.t); % 34 percent for M77

M77_saccade = importdata('YutaTest77c_REMsWAKE.mat');
%% Get pupil position data, AWAKE
% M77_pupil = importdata('YutaTest77c_eye_pupil_positions_converted.mat'); % pupil positions during wakefullness. Don't need this (yet), just saccades
% M77_saccade.CommonStart are the saccade start times, in seconds. ***I think that all saccades (for both eyes) are concatenated
% M77_saccade.ampX{1}: horizontal amplitude of saccades in RIGHT eye, in degree.
% M77_saccade.ampX{2}: horizontal amplitude of rapid saccades in LEFT eye, in degree.

M77_saccade.t = M77_saccade.CommonStart;
M77_num_saccades = length(M77_saccade.t); % 9065 for M77
M77_saccade.right = M77_saccade.ampX{1}; % right saccade amplitudes
M77_saccade.left = M77_saccade.ampX{2};  % left saccade amplitudes

cfg_in = [];
cfg_in.window = [-.2 .2];
doPlot = 1;
startCell = 1; endCell = length(M77.t);
if doPlot
    for iCell = startCell:endCell
        disp(num2str(iCell))
        clf
        myCell = SelectTS([], M77, iCell);
        myCell.type = 'ts';
        tic
        [outputS, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, M77_saccade.t, 'doPlot', doPlot);
        toc
        %     [peth_out, all_trials] = TSDpeth_fast(cfg_in, myCell, M77_saccade.t); % I don't think that I can use tsdPETH to do a spikePETH
        title(num2str(iCell))
        %     disp('press any key to continue')
        pause
    end
end
%% Isolate times for each neuron when the mouse is near the PFD (preferred firing direction) 











% ------------------------------------------------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------------------------------------------------------------------
%% M79
%% Import the spike train data
disp('importing data')
M79_TCs = importdata('YutaTest79c_HeadDirection_OpenField.mat'); % data strcuture that contains Maps, kappa_ppln, pval_ppln, z_ppln.
% Maps contains the tuning curves. pval_ppln is presumably a p value for some significance test for head direction selectivity.
% Not sure what kappa is. z_ppln might be z-score.
numCellsM79 = length(M79_TCs.Maps);

spiketrain_1_M79 = importdata('YutaTest79c.res.1'); % Shank 1 from mouse 77. I think these spike trains are split in two because of their size.
spiketrain_1_M79 = spiketrain_1_M79./20000; % divide by the sampling rate (20kHz) to get seconds
spiketrain_2_M79 = importdata('YutaTest79c.res.2'); % Shank 2 from mouse 77.
spiketrain_2_M79 = spiketrain_2_M79./20000; % divide by the sampling rate (20kHz) to get seconds


clu_1_M79 = importdata('YutaTest79c.clu.1');  % 'clu_1' here means the clusters from shank 1. 'clu_2' means the clusters from shank 2.
M79shank1_num = clu_1_M79(1); % this includes zeros and ones
clu_1_M79toUse = clu_1_M79(2:end);
clu_2_M79 = importdata('YutaTest79c.clu.2');
M79shank2_num = clu_2_M79(1); % this includes zeros and ones
clu_2_M79toUse = clu_2_M79(2:end);

assert(length(spiketrain_1_M79) == length(clu_1_M79toUse))
assert(length(spiketrain_2_M79) == length(clu_2_M79toUse))
%% Isolate the individual clusters (i.e., the single units)
clu_1_M79_IDs = unique(clu_1_M79toUse);
temp1 = ismember(clu_1_M79_IDs, 0:1); % 0 = artifact and 1 = noise. Remove these from consideration.
num_clu_1_M79 = sum(temp1==0); % how many single units (clusters) are there?
clu_1_M79_IDsToUse = clu_1_M79_IDs(clu_1_M79_IDs ~= 0 & clu_1_M79_IDs ~= 1); % Use clusters that are not zero or one.

clu_2_M79_IDs = unique(clu_2_M79toUse);
temp2 = ismember(clu_2_M79_IDs, 0:1); % 0 = artifact and 1 = noise. Remove these from consideration.
num_clu_2_M79 = sum(temp2==0); % how many single units (clusters) are there?
clu_2_M79_IDsToUse = clu_2_M79_IDs(clu_2_M79_IDs ~= 0 & clu_2_M79_IDs ~= 1); % Use clusters that are not zero or one.

totalCellsM79 = num_clu_1_M79 + num_clu_2_M79;
if numCellsM79 ~= totalCellsM79
    warning('M77 neurons do not add up')
end

%% Create the spike train structure
clear M79
counterToUse = 0;
M79.type = 'ts';
for iC = 1: totalCellsM79  % combine spike trains from each file into one structure
    if iC <= num_clu_1_M79
        M79.t{iC} = spiketrain_1_M79(clu_1_M79toUse == clu_1_M79_IDsToUse(iC));
    elseif iC > num_clu_1_M79 && iC <= totalCellsM79
        counterToUse = counterToUse + 1;
        M79.t{iC} = spiketrain_2_M79(clu_2_M79toUse == clu_2_M79_IDsToUse(counterToUse));
    else
        warning('problem with neuron count')
    end
    M79.label{iC} = 'M79';
end


%% Get heading and AHV
M79_heading = importdata('YutaTest79c_OpenField_HeadDirection.mat'); % contains Neck&NoseOmitIdx [value of 1 = omit], rho(distance between nose to neck),
% theta (heading, in radians), t (timestamp, IN SECONDS)
M79tracking_dur = (M79_heading.t(end) - M79_heading.t(1))/60; % How many minutes long the tracking was in the open field. This is 98 minutes for M77.
M79_Q = unwrap(M79_heading.theta);
M79_Q_deg = M79_Q*180/pi;
window = 0.1;
postsmoothing = .05;
M79_AHV = dxdt(M79_heading.t, M79_Q_deg, 'window', window, 'postsmoothing', postsmoothing);
M79_AHVtsd = tsd(M79_heading.t, M79_AHV);
M79_heading_sampling_rate = 1/median(diff(M79_heading.t));  % should be 50Hz
% Isolate low AHV times
AHVthresh = 10; % cm/sec
M79_low_AHV = abs(M79_AHV) < AHVthresh;
M79_low_AHV_samples = sum(M79_low_AHV);
M79_low_AHV_seconds = M79_low_AHV_samples/M79_heading_sampling_rate;
M79_low_AHV_minutes = M79_low_AHV_seconds/60;
fraction_M79_low_AHV = M79_low_AHV_samples / length(M79_heading.t); % 34 percent for M77

M79_saccade = importdata('YutaTest79c_REMsWAKE.mat');
%% Get pupil position data, AWAKE
% M77_pupil = importdata('YutaTest77c_eye_pupil_positions_converted.mat'); % pupil positions during wakefullness. Don't need this (yet), just saccades
% M77_saccade.CommonStart are the saccade start times, in seconds. ***I think that all saccades (for both eyes) are concatenated
% M77_saccade.ampX{1}: horizontal amplitude of saccades in RIGHT eye, in degree.
% M77_saccade.ampX{2}: horizontal amplitude of rapid saccades in LEFT eye, in degree.

M79_saccade.t = M79_saccade.CommonStart;
M79_num_saccades = length(M79_saccade.t); % 9065 for M77
M79_saccade.right = M79_saccade.ampX{1}; % right saccade amplitudes
M79_saccade.left = M79_saccade.ampX{2};  % left saccade amplitudes

cfg_in = [];
cfg_in.window = [-.2 .2];
doPlot = 1;
startCell = 1; endCell = length(M79.t);
if doPlot
    for iCell = startCell:endCell
        disp(num2str(iCell))
        clf
        myCell = SelectTS([], M79, iCell);
        myCell.type = 'ts';
        tic
        [outputS, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, M79_saccade.t, 'doPlot', doPlot);
        toc
        %     [peth_out, all_trials] = TSDpeth_fast(cfg_in, myCell, M77_saccade.t); % I don't think that I can use tsdPETH to do a spikePETH
        title(num2str(iCell))
        %     disp('press any key to continue')
        pause
    end
end








%% extra
% %% Determine how many neurons are HD cells, based on the p-values provided
% sig77 = pval_ppln_77 < .05; num_sig77 = sum(sig77); fract_sig_77 = num_sig77/length(Maps77);
% sig79 = pval_ppln_79 < .05; num_sig79 = sum(sig79); fract_sig_79 = num_sig79/length(Maps79);
% sig83 = pval_ppln_83 < .05; num_sig83 = sum(sig83); fract_sig_83 = num_sig83/length(Maps83);
% sig85 = pval_ppln_85 < .05; num_sig85 = sum(sig85); fract_sig_85 = num_sig85/length(Maps85);
% sig96 = pval_ppln_96 < .05; num_sig96 = sum(sig96); fract_sig_96 = num_sig96/length(Maps96);
% sig98 = pval_ppln_98 < .05; num_sig98 = sum(sig98); fract_sig_98 = num_sig98/length(Maps98);
%
% open_field_sig_fraction = [fract_sig_77 fract_sig_79 fract_sig_83 fract_sig_85 fract_sig_96 fract_sig_98];