% 2024-10-24. JJS. Extracts Senzai data for one mouse.
TCbinCenters = 5:10:355; 
cd('J:\senzai dataset\data files')
%% M77
%% Import the spike train data
disp('importing data')
M77_TCs = importdata('YutaTest77c_HeadDirection_OpenField.mat'); % data strcuture that contains Maps, kappa_ppln, pval_ppln, z_ppln.
% Maps contains the tuning curves. pval_ppln is presumably a p value for some significance test for head direction selectivity.
% Not sure what kappa is. z_ppln might be z-score.
numCellsM77 = length(M77_TCs.Maps);

spiketrain_1_M77 = importdata('YutaTest77c.res.1'); % Shank 1 from mouse 77. 
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
temp1 = ismember(clu_1_M77_IDs, 0:1); % 0 = artifact and 1 = noise. Remove these from consideration.git 
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

% **************************************************************************************************************************************************************
% **************************************************************************************************************************************************************
%% Get heading and AHV
M77_heading = importdata('YutaTest77c_OpenField_HeadDirection.mat'); % contains Neck&NoseOmitIdx [value of 1 = omit], rho(distance between nose to neck),
% theta (heading, in radians), t (timestamp, IN SECONDS)
M77_heading_degrees = M77_heading.theta*180/pi;
M77_heading_degreesTSD = tsd(M77_heading.t, M77_heading_degrees);

M77tracking_dur = (M77_heading.t(end) - M77_heading.t(1))/60; % How many minutes long the tracking was in the open field. This is 98 minutes for M77.
M77_Q = unwrap(M77_heading.theta);
M77_Q_deg = M77_Q*180/pi;
window = 0.1;
postsmoothing = .05;
M77_AHV = dxdt(M77_heading.t, M77_Q_deg, 'window', window, 'postsmoothing', postsmoothing);
M77_AHVtsd = tsd(M77_heading.t, M77_AHV);
M77_heading_sampling_rate = 1/median(diff(M77_heading.t));  % should be 50Hz
% **************************************************************************************************************************************************************
% **************************************************************************************************************************************************************
%% Isolate low AHV times
% Get Indices
AHVthresh = 10; % cm/sec
M77_low_AHV = abs(M77_AHV) < AHVthresh;
M77_low_AHV_samples = sum(M77_low_AHV);
M77_low_AHV_seconds = M77_low_AHV_samples/M77_heading_sampling_rate;
M77_low_AHV_minutes = M77_low_AHV_seconds/60;
fraction_M77_low_AHV = M77_low_AHV_samples / length(M77_heading.t); % 34 percent for M77
% Get Start/Stop times (for later restrict)

M77_low_AHV_diff = horzcat(NaN, diff(M77_low_AHV)); 
% High to Low
% ---------------------------------------
%                                       -
%                                       -
%                                       -
%                                       -
%                                       ------------------------------------------------
low_AHV_ones = find(M77_low_AHV_diff == 1); sum_low_AHV_ones = length(low_AHV_ones); % find the transition points from high(er) AHV to low AHV 
low_AHV_tstart = M77_AHVtsd.tvec(low_AHV_ones); low_AHV_tstart = low_AHV_tstart';
% Low to High 
%                                       ------------------------------------------------
%                                       -
%                                       -
%                                       -
%                                       -
% ---------------------------------------
low_AHV_minus_ones = find(M77_low_AHV_diff == -1); sum_low_AHV_minus_ones = length(low_AHV_minus_ones);  % find transitions from low AHV back to high AHV
low_AHV_tend = M77_AHVtsd.tvec(low_AHV_minus_ones); low_AHV_tend = low_AHV_tend';

% Account for the first value being low or high. *** Need to account for all possiblities *** 

if low_AHV_tstart(1) > low_AHV_tend(1)   % In other words, if the first transition is from high to low, then the session started out LOW. Make tstart(1) = the first timestamp
    low_AHV_tstart = [M77_AHVtsd.tvec(1) low_AHV_tstart];
end
if low_AHV_tstart(end) > low_AHV_tend(end) % In other words, if the last transition is from high to low, then the session ended LOW. Make tend(end) = last timestamp
    low_AHV_tend = [low_AHV_tend M77_AHVtsd.tvec(end)];
end
assert(length(low_AHV_tstart) == length(low_AHV_tstart))

%% Load the Saccade data
M77_saccade = importdata('YutaTest77c_REMsWAKE.mat');
% M77_pupil = importdata('YutaTest77c_eye_pupil_positions_converted.mat'); % pupil positions during wakefullness. Don't need this (yet), just saccades
% M77_saccade.CommonStart are the saccade start times, in seconds. ***I think that all saccades (for both eyes) are concatenated
% M77_saccade.ampX{1}: horizontal amplitude of saccades in RIGHT eye, in degree.
% M77_saccade.ampX{2}: horizontal amplitude of rapid saccades in LEFT eye, in degree.

M77_saccade.t = M77_saccade.CommonStart;
M77_num_saccades = length(M77_saccade.t); % 9065 for M77
M77_saccade.right = M77_saccade.ampX{1}; % right saccade amplitudes
M77_saccade.left = M77_saccade.ampX{2};  % left saccade amplitudes

% cfg_in = [];
% cfg_in.window = [-.2 .2];
% doPlot = 1;
% startCell = 1; endCell = length(M77.t);
% % Awake saccade peths, all saccades
% if doPlot
%     for iCell = startCell:endCell
%         disp(num2str(iCell))
%         clf
%         myCell = SelectTS([], M77, iCell);
%         myCell.type = 'ts';
%         tic
%         [outputS, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, M77_saccade.t, 'doPlot', doPlot);
%         toc
%         %     [peth_out, all_trials] = TSDpeth_fast(cfg_in, myCell, M77_saccade.t); % I don't think that I can use tsdPETH to do a spikePETH
%         title(num2str(iCell))
%         %     disp('press any key to continue')
%         pause
%     end
% end

%% Generate saccade-triggered AHV PETH                  ...this isn't very interesting will all the data included b/c its just random and noisy. 
cfg_in.mode = 'interp';
cfg_in.window = [-.5 .5];
cfg_in.dt = .02;
if doPlot
    [peth_out, all_trials] = TSDpeth_fast(cfg_in, M77_AHVtsd, M77_saccade.t); % I don't think that I can use tsdPETH to do a spikePETH
end
a=colorbar;
a.Label.String = 'AHV (deg/s)';

%% Choose neurons with minimum FR and Strong HD tuning   [M77_neuronsToUse]
% Find the PFD (peak) for each neuron that has significant HD tuning
HD_z_thresh = 50;  % z-score threshold. ***Not sure how Yuta calculates the z-score. Check on this.
% Peak firing rate criterion
for iC = 1:length(M77_TCs.Maps)
    M77_maxFR(iC) = max(M77_TCs.Maps{iC,1}.rate);
    M77_maxFRsmooth(iC) = max(smoothdata(M77_TCs.Maps{iC,1}.rate));
end
FRthresh = 10; % in Hz. The tuning curve peak must be at or above this threshold to include the neuron. 
M77_peakFRtoUse = M77_maxFR > FRthresh;
M77_peakFRtoUseIDs = find(M77_peakFRtoUse); M77_peakFRtoUseIDs = M77_peakFRtoUseIDs';
sum_M77_peakFRtoUse = sum(M77_peakFRtoUse); fraction_M77_peakFRtoUse = sum_M77_peakFRtoUse/length(M77_TCs.Maps);

M77_HDC_logical = M77_TCs.z_ppln > HD_z_thresh;
M77_HDC_IDs = find(M77_HDC_logical);
M77_sum_HDCs = sum(M77_HDC_logical);

M77_neuronsToUse = intersect(M77_peakFRtoUseIDs, M77_HDC_IDs);
sum_M77_HDC_minFR = length(M77_neuronsToUse);

clf; hold on; set(gca, 'FontSize', 24)
doNorm = 0;
if doPlot
    for iC = 1 : sum_M77_HDC_minFR
        IDtoUse(iC) = M77_neuronsToUse(iC);
        disp(num2str(IDtoUse(iC)))
        normIC(iC,:) = (M77_TCs.Maps{iC,1}.rate - min(M77_TCs.Maps{iC,1}.rate)) ./ (max(M77_TCs.Maps{iC,1}.rate) - min(M77_TCs.Maps{iC,1}.rate));
        M77_TCs.Maps{iC,1}.rateSmoothed = smoothdata(M77_TCs.Maps{iC,1}.rate);
        if doNorm
            plot(TCbinCenters, smoothdata(normIC(iC,:)));
            ylabel('normalized FR (Hz)')
        else
            plot(TCbinCenters, smoothdata(M77_TCs.Maps{IDtoUse(iC),1}.rate));
%             line([halfwidth(iC,1) halfwidth(iC,1)], [0 M77_TCs_smoothed(iC, index1(iC))], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
%             line([halfwidth(iC,2) halfwidth(iC,2)], [0 M77_TCs_smoothed(iC, index2(iC))], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
            ylabel('FR (Hz)')
        end
        title(strcat('M77 cell num', num2str(IDtoUse)))
        xlabel('HD (deg)')
        pause
    end
end
c = axis;
axis([0 365 c(3) c(4)]);
%% AWAKE Saccade PETHs, restricted to low AHV
% low_AHV_tstart and low_AHV_tend are the    start and end times for low AHV intervals 
for iC = 
    
    
    
    
    
    
    
    
    
    
    

%% Restrict to PFD times
% Isolate the PFD peak
% Estimate PFD range  [half maximum]
for iC = 1 : sum_M77_HDC_minFR
    M77_TCs_smoothed(iC,:) = smoothdata(M77_TCs.Maps{IDtoUse(iC),1}.rate);    
%     [peakFR(iC), peakIndex(iC)] = max(M77_maxFRsmooth(IDtoUse(iC)));
    [minFR(iC), minIndex(iC)] = min(M77_TCs_smoothed(iC,:)); minIndexDeg(iC) = TCbinCenters(minIndex(iC));
    [maxFR(iC), maxIndex(iC)] = max(M77_TCs_smoothed(iC,:)); maxIndexDeg(iC) = TCbinCenters(maxIndex(iC));
    halfMax(iC) = (minFR(iC) + maxFR(iC))/2;
    
    % *** next 4 lines are glitchy and don't work right if the first bin is the lowest or highest value ***
    % Find where the data first drops below half the max.
    index1(iC) = find(M77_TCs_smoothed(iC,:) >= halfMax(iC), 1, 'first'); 
    halfwidth(iC,1) = TCbinCenters(index1(iC));
    % Find where the data last rises above half the max.
    index2(iC) = find(M77_TCs_smoothed(iC,:) >= halfMax(iC), 1, 'last');  
    halfwidth(iC,2) = TCbinCenters(index2(iC));

    above(iC,:) = M77_TCs_smoothed(iC,:) > halfMax(iC);
    range{iC} = find(above(iC,:));
    
    FWHM{iC} = TCbinCenters(range{iC});
end
% Plot the chosen TCs with the min, max, and half-max lines 
for iC = 1 : sum_M77_HDC_minFR
    clf; hold on
    disp(num2str(IDtoUse(iC)))
    plot(TCbinCenters, M77_TCs_smoothed(iC,:));
    line([minIndexDeg(iC) minIndexDeg(iC)], [0 M77_TCs_smoothed(iC, minIndex(iC))], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1) % min line (often not visible)
    line([maxIndexDeg(iC) maxIndexDeg(iC)], [0 M77_TCs_smoothed(iC, maxIndex(iC))], 'LineStyle', '--', 'Color', 'g', 'LineWidth', 1) % max line (green)
    
    line([halfwidth(iC,1) halfwidth(iC,1)], [0 M77_TCs_smoothed(iC, index1(iC))], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
    line([halfwidth(iC,2) halfwidth(iC,2)], [0 M77_TCs_smoothed(iC, index2(iC))], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
    
    ylabel('FR (Hz)')
    title(num2str(IDtoUse(iC)))
    pause
end
a=colorbar;
a.Label.String = 'FR';
%% 






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