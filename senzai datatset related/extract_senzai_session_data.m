% 2024-10-24. JJS. Extracts Senzai data for one mouse.
TCbinCenters = -175:10:175;
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
M77_heading_degreesTSD = tsd(M77_heading.t, M77_heading_degrees');

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

%% RESTRICT THE DATA TO LOW AHV TIMES
M77_heading_degreesTSD_lowAHV = restrict(M77_AHVtsd, low_AHV_tstart, low_AHV_tend);

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

% These threshold values are specific to M77
clf; plot(M77_saccade.right, M77_saccade.left, '.'); hold on
posRightLarge = M77_saccade.right > 5;
posLeftLarge = M77_saccade.left > 12;
posLeftRightLarge = posRightLarge & posLeftLarge; NUMposLeftRightLarge = sum(posLeftRightLarge);
plot(M77_saccade.right(posLeftRightLarge), M77_saccade.left(posLeftRightLarge), 'r.')

negRightLarge = M77_saccade.right < -11;
negLeftLarge = M77_saccade.left < -6;
negLeftRightLarge = negRightLarge & negLeftLarge; NUMnegLeftRightLarge = sum(negLeftRightLarge);
bothLarge = posLeftRightLarge | negLeftRightLarge; sum_bothLarge = sum(bothLarge); disp(sum_bothLarge);
LargeSaccades = M77_saccade.t(bothLarge);

% ***************************** Double check the terminology. Does a positive amplitdue saccade = CW? ***************************
CCW_large_saccades = M77_saccade.t(negLeftRightLarge);
CW_large_saccades = M77_saccade.t(posLeftRightLarge);

plot(M77_saccade.right(negLeftRightLarge), M77_saccade.left(negLeftRightLarge), 'm.')
text(-30, 30, strcat('num Pos =', num2str(NUMposLeftRightLarge)));
text(-30, 25, strcat('num Neg =', num2str(NUMnegLeftRightLarge)));
text(-30, 20, strcat('total =', num2str(NUMposLeftRightLarge + NUMnegLeftRightLarge)));
set(gca, 'FontSize', 26)
xlabel('right eye amplitude')
ylabel('left eye amplitude')

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
% doPlot = 1;
% cfg_in.mode = 'interp';
% cfg_in.window = [-.5 .5];
% cfg_in.dt = .02;
% if doPlot
%     [peth_out, all_trials] = TSDpeth_fast(cfg_in, M77_AHVtsd, M77_saccade.t); % I don't think that I can use tsdPETH to do a spikePETH
% end
% a=colorbar;
% a.Label.String = 'AHV (deg/s)';

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
% axis([0 365 c(3) c(4)]);
%% AWAKE Saccade PETHs, restricted to low AHV.  [M77_heading_degreesTSD_lowAHV]

low_ahv_Interval = iv(low_AHV_tstart, low_AHV_tend);
large_saccades_ts = ts({LargeSaccades});
large_amplitude_low_ahv_saccades = restrict(large_saccades_ts, low_ahv_Interval);
cfg_in = [];
cfg_in.window = [-.2 .2];
doPlot = 1;
startCell = 1; endCell = length(M77.t);
% Awake saccade peths, all saccades
if doPlot
    for iCell = startCell:endCell
        disp(num2str(iCell))
        clf
        myCell = SelectTS([], M77, iCell);
        myCell.type = 'ts';
        tic
        [outputS, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, large_amplitude_low_ahv_saccades.t{1}, 'doPlot', doPlot);
        toc
        title(num2str(iCell))
        pause
    end
end

%% Restrict to PFD times
% Isolate the PFD peak
% Estimate PFD range  [half maximum (or whatever fraction)]
fractiontouse = 0.2;   % lower the value here, the more data is included (more of the tuning curve)
rangeNaNs = NaN(length(M77_neuronsToUse),36);
for iC = 1 : sum_M77_HDC_minFR
    M77_TCs_smoothed(iC,:) = smoothdata(M77_TCs.Maps{IDtoUse(iC),1}.rate);
    %     [peakFR(iC), peakIndex(iC)] = max(M77_maxFRsmooth(IDtoUse(iC)));
    [minFR(iC), minIndex(iC)] = min(M77_TCs_smoothed(iC,:)); minIndexDeg(iC) = TCbinCenters(minIndex(iC));
    [maxFR(iC), maxIndex(iC)] = max(M77_TCs_smoothed(iC,:)); maxIndexDeg(iC) = TCbinCenters(maxIndex(iC));
    %     halfMax(iC) = (minFR(iC) + maxFR(iC))/2;
    halfMax(iC) = minFR(iC) + ( maxFR(iC) - minFR(iC))*fractiontouse;
    
    % *** next 4 lines are glitchy and don't work right if the first bin is the lowest or highest value ***
    % Find where the data first drops below half the max.
    index1(iC) = find(M77_TCs_smoothed(iC,:) >= halfMax(iC), 1, 'first');
    halfwidth(iC,1) = TCbinCenters(index1(iC));
    % Find where the data last rises above half the max.
    index2(iC) = find(M77_TCs_smoothed(iC,:) >= halfMax(iC), 1, 'last');
    halfwidth(iC,2) = TCbinCenters(index2(iC));
    
    above(iC,:) = M77_TCs_smoothed(iC,:) > halfMax(iC);
    range{iC} = find(above(iC,:));
    rangeNaNs(iC,range{iC}) = TCbinCenters(range{iC});
    
    FWHM{iC} = TCbinCenters(range{iC});   % **** this gives the wrong answer when the TC peak is at the edge(s)
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
%% Break up Saccades into CW or CCW of PFD, and in or out of field    [CW_large_saccades & CCW_large_saccades]
% define the intervals when the mouse is
%       (1) in the FWHM of the tuning curve
%       (2) to the left (or right) of the peak
%           and (3) makes a saccade TOWARD (or AWAY) from the peak
% LEFT -> CW = toward. LEFT -> CCW = away. RIGHT -> CCW = toward. RIGHT -> = AWAY.

% Example:  FWHM{1,2} = [175   185   195   205   215   225   235   245   255   265]; Need to select headings >= 170 & <= 270

% M77_neuronsToUse = intersect(M77_peakFRtoUseIDs, M77_HDC_IDs);
% IDtoUse(iC) = M77_neuronsToUse(iC);
% M77_Q = unwrap(M77_heading.theta);
% M77_Q_deg = M77_Q*180/pi;
% M77_heading_degrees
clear rangetouse
%  from 0 degrees to first non-NaN bin, OR from last non-NaN bin to 360 degrees.
for iNeuron = 1: length(M77_neuronsToUse)
    fwhm = FWHM{1, iNeuron};
    rangetouse(iNeuron,1) = min(fwhm) - 5; % remember that bin centers are at intervals of 5,15,25, etc. BinEdges are 0,10,20, etc.
    rangetouse(iNeuron,2) = max(fwhm) + 5;
    
    if rangetouse(iNeuron,1) == -180 && rangetouse(iNeuron,2) == 180
        wrapAroundCell(iNeuron) = 1;
        %         rangetouse(iNeuron,:) = NaN;
        rangetouse(iNeuron,1) = NaN;
        rangetouse(iNeuron,2) = NaN;
    else
        wrapAroundCell(iNeuron) = 0;
    end
end
wrapAroundCell = wrapAroundCell';

%% For each Neuron, restrict to when mouse is in PFD
% ******************************************************** Remember that heading is between -180 and +180 here, not 0 - 360 ****************************************
samples = length(M77_heading_degreesTSD.tvec);
WRAPrangetouse = NaN(length(M77_neuronsToUse),4);
for iNeuron = 1: length(M77_neuronsToUse)
    if wrapAroundCell(iNeuron) == 0         % 'easy' cells
        M77_headings_to_use{iNeuron} =  M77_heading_degreesTSD.data > rangetouse(iNeuron,1) & M77_heading_degreesTSD.data < rangetouse(iNeuron,2);
    elseif wrapAroundCell(iNeuron) == 1     % wraparound cells
        temp = isnan(rangeNaNs(iNeuron,:));
        a = find(temp); % find the first empty bin
        if a(1) == 1
            b = 1;
        else
            b = a(1) - 1;
        end
        if a(end) == 36
            b = 36;
        else
            c = a(end) + 1;
        end
        WRAPrangetouse(iNeuron,1) = -180;
        WRAPrangetouse(iNeuron,2) = TCbinCenters(b) + 5; % remember that bin centers are at intervals of 5,15,25, etc. BinEdges are 0,10,20, etc.
        WRAPrangetouse(iNeuron,3) = TCbinCenters(c) - 5;
        WRAPrangetouse(iNeuron,4) = 180;
        
        M77_headings_to_use{iNeuron} = M77_heading_degreesTSD.data < WRAPrangetouse(iNeuron,2) | M77_heading_degreesTSD.data > WRAPrangetouse(iNeuron,3);
        %  **** Need to turn this into tStart and tEnd times
    else
        error('problem')
    end
    fraction_left(iNeuron) = sum(M77_headings_to_use{iNeuron})/samples;
end
fraction_left = fraction_left'; if doPlot; clf; plot(fraction_left); title('fraction of the session left after restricting to PFD'); end
%% Get the Peaks
for iNeuron = 1:length(M77_neuronsToUse)
    [peakValue(iNeuron), peakIndex(iNeuron)] = max(smoothdata(M77_TCs.Maps{M77_neuronsToUse(iNeuron)}.rate));  % peakIndex = which bin in TCbinCenters is the peak firing rate for that neuron
end
rangetouse_with_peak(:,1) = rangetouse(:,1);
rangetouse_with_peak(:,2) = TCbinCenters(peakIndex)';
rangetouse_with_peak(:,3) = rangetouse(:,2);
tempNotNaN = ~isnan(rangetouse(:,1));
tempNaN = isnan(rangetouse(:,1));
HDnums = 1:length(M77_HDC_IDs);
nonWrapCells = HDnums(tempNotNaN);
WrapCells = HDnums(tempNaN);
% rangetouse_with_peak = rangetouse_with_peak(nonWrapCells,:);  % 12 x 3


%% Find samples when the mouse is CW or CCW to the peak.
% **********check this****************    CCW = (+), CW = (-)
clear tCW; clear tCCW; clear entering_CW_zone; clear exiting_CW_zone; clear start_stop_table_CW; clear ssInterval_CW; clear entering_CW_zone_timestamps; clear exiting_CW_zone_timestamps; clear CW_large_saccades_ts; clear CW_CW
tCW = cell(1,length(M77_neuronsToUse));   % Logical where NaNs = not in PFD 
tCCW = cell(1,length(M77_neuronsToUse));
numHeadingSamples = length(M77_heading_degreesTSD.tvec);
heading_sampling_rate = 1/median(diff(M77_heading_degreesTSD.tvec));  % this should be 50 Hz. Each sample = 20ms.  
% *** Need to add a restriction for cases where the interval is too short, like, less than 200ms?
for iNeuron = 1:length(M77_neuronsToUse)
    if ismember(iNeuron, WrapCells)
        tCW{iNeuron} = {};
        tCCW{iNeuron} = {};
        ssInterval_CW{iNeuron} = {};
        entering_CW_zone_timestamps{iNeuron} = {};
        exiting_CW_zone_timestamps{iNeuron} = {}; 
        CW_CW{iNeuron}= {};
        CCW_CCW{iNeuron} = {}; 
    else
        tCW{iNeuron} = M77_heading_iNeuron{iNeuron}.data > rangetouse_with_peak(iNeuron,1) & M77_heading_iNeuron{iNeuron}.data < rangetouse_with_peak(iNeuron,3);
        tCCW{iNeuron} = M77_heading_iNeuron{iNeuron}.data < rangetouse_with_peak(iNeuron,3) & M77_heading_iNeuron{iNeuron}.data > rangetouse_with_peak(iNeuron,2);
        
        diff_tCW{iNeuron} = diff(tCW{iNeuron}); diff_tCW{iNeuron} = horzcat(NaN, diff_tCW{iNeuron}); 
        diff_tCCW{iNeuron} = diff(tCCW{iNeuron}); diff_tCCW{iNeuron} = horzcat(NaN, diff_tCCW{iNeuron}); 
        
        entering_CW_zone{iNeuron} = find(diff_tCW{iNeuron} == 1); 
        exiting_CW_zone{iNeuron} = find(diff_tCW{iNeuron} == -1); 
        
        start_stop_table_CW(iNeuron,1) = length(entering_CW_zone{iNeuron}); 
        start_stop_table_CW(iNeuron,2) = length(exiting_CW_zone{iNeuron}); 
        
        if start_stop_table_CW(iNeuron,1) > start_stop_table_CW(iNeuron,2)   % in this case, session ended in CW zone. 1 'extra' start time. Make last stop time t(end)
            exiting_CW_zone{iNeuron} = horzcat(exiting_CW_zone{iNeuron}, numHeadingSamples); % numHeadingSamples is also the index of the last sample
        end
        if start_stop_table_CW(iNeuron,2) > start_stop_table_CW(iNeuron,1)   % in this case, session started in CW zone. 1 'extra' end time. Make first stop time t(1)
            entering_CW_zone{iNeuron} = horzcat(entering_CW_zone{iNeuron}, 1); % numHeadingSamples is also the index of the last sample
        end  
        
        start_stop_table_CW(iNeuron,1) = length(entering_CW_zone{iNeuron}); 
        start_stop_table_CW(iNeuron,2) = length(exiting_CW_zone{iNeuron}); 
        assert(start_stop_table_CW(iNeuron,1) == start_stop_table_CW(iNeuron,2))
        ssInterval_CW{iNeuron} = exiting_CW_zone{iNeuron} - entering_CW_zone{iNeuron}; 
        
        entering_CW_zone_timestamps{iNeuron} = M77_heading_degreesTSD.tvec(entering_CW_zone{iNeuron});
        exiting_CW_zone_timestamps{iNeuron} = M77_heading_degreesTSD.tvec(exiting_CW_zone{iNeuron});
        
        CW_large_saccades_ts = ts({CW_large_saccades});  % needs to be a ts structure in order to use restrict
        CCW_large_saccades_ts = ts({CCW_large_saccades});
        
        CW_CW{iNeuron} = restrict(CW_large_saccades_ts, entering_CW_zone_timestamps{iNeuron}, exiting_CW_zone_timestamps{iNeuron});  % CW zone, CW saccade. TOWARD
        CW_CCW{iNeuron} = restrict(CCW_large_saccades_ts, entering_CW_zone_timestamps{iNeuron}, exiting_CW_zone_timestamps{iNeuron}); % CW zone, CCW saccade. AWAY

        CW_saccade_table(iNeuron,1) = length(CW_CW{iNeuron}.t{1});   % TOWARD 
        CW_saccade_table(iNeuron,2) = length(CW_CCW{iNeuron}.t{1});  % AWAY 
    end                                                                                  % stands for start-stop interval
end
%% Plot the saccade PETHs

cfg_in = [];
cfg_in.window = [-.2 .2];
doPlot = 1;
startCell = 1; endCell = length(M77_neuronsToUse);
% Awake saccade peths, all saccades
if doPlot
    for iCell = 1:length(M77_neuronsToUse)   % M77_neuronsToUse(iCell) 
        disp(num2str(iCell))
        clf
        myCell = SelectTS([], M77, M77_neuronsToUse(iCell));
        myCell.type = 'ts';
        tic
        [outputS, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, CW_CW{iNeuron}.t{1}, 'doPlot', doPlot);
        title('CW zone, CW saccade - TOWARD') 
        pause
        [outputS, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, CW_CCW{iNeuron}.t{1}, 'doPlot', doPlot);
        title('CW zone, CCW saccade - AWAY')
        toc
%         title(num2str(iCell))
        pause
    end
end

















%% Plot the CW and CCW intervals to make sure that they look right
for iNeuron = 1:length(M77_neuronsToUse)    % plot the Heading data with overlay of the restricted part for each neuron's PFD
    clf
    if ~ismember(iNeuron, WrapCells)
        plot(M77_heading_degreesTSD.tvec, M77_heading_degreesTSD.data, '.'); hold on
        %     plot(M77_heading_degreesTSD.tvec(M77_headings_to_use{iNeuron}), M77_heading_degreesTSD.data(M77_headings_to_use{iNeuron}), '.');
        plot(M77_heading_iNeuron{iNeuron}.tvec(tCW{iNeuron}), M77_heading_iNeuron{iNeuron}.data(tCW{iNeuron}), 'r.')    % CW 
        plot(M77_heading_iNeuron{iNeuron}.tvec(tCCW{iNeuron}), M77_heading_iNeuron{iNeuron}.data(tCCW{iNeuron}), 'g.')  % CCW 
        title(num2str(iNeuron))   % how much of the TC is included = 1 - fractiontouse. If fractiontouse = 0.2, then data is restricted to that 80% of TC.
        % If full tuning curve with is 180 deg, then a 0.2 value would mean that we are keeping .8(360-180/360)  40% of the overall session data.
        xlabel('Heading')
        ylabel('Heading')
        a = xlim;
        xrange = a(2) - a(1);
        xUnit = xrange/36;
        binRange = 360;
        newBinCenters = a(1)+0.5*(xUnit):xUnit:a(2)-0.5*(xUnit);
        
        yyaxis right
        plot(newBinCenters, smoothdata(M77_TCs.Maps{M77_neuronsToUse(iNeuron)}.rate), 'LineWidth', 15)
        xticks(newBinCenters);
        xticklabelstouse = {'-175','-165','-155','-145','-135','-125','-115','-105','-95','-85','-75','-65','-55','-45','-35','-25','-15','-5','5','15','25','35','45','55','65','75','85','95','105','115','125','135','145','155','165','175'};
        xticklabels(xticklabelstouse)
        set(gca, 'FontSize', 22)
        line([newBinCenters(index1(iNeuron)) newBinCenters(index1(iNeuron))], [0 halfMax(iNeuron)], 'LineWidth', 15, 'Color', 'Red')
        line([newBinCenters(index2(iNeuron)) newBinCenters(index2(iNeuron))], [0 halfMax(iNeuron)], 'LineWidth', 15, 'Color', 'Green')
        line([newBinCenters(maxIndex(iNeuron)) newBinCenters(maxIndex(iNeuron))], [0 maxFR(iNeuron)], 'LineWidth', 15)
        saveas(gcf, strcat(num2str(iNeuron),'.png'));
        y = ylim;
        text(newBinCenters(6), .87*(y(2)), 'CW', 'FontSize', 55, 'Color', 'Red');
        text(newBinCenters(28), .87*(y(2)), 'CCW', 'FontSize', 55, 'Color', 'Green');
        disp('saving fig')
        saveas(gcf,strcat(num2str(iNeuron), '.png'))        
        pause
    end
end

%% For each neuron, get the START/STOP times for Headings within the PFD epochs
% Additionally, separate out by CW of PFD peak and CCW of PFD peak

for iNeuron = 1:length(M77_neuronsToUse)
    M77_heading_iNeuron{iNeuron} = M77_heading_degreesTSD;
    NaNindex = M77_headings_to_use{iNeuron} == 0;
    M77_heading_iNeuron{iNeuron}.data(NaNindex) = NaN;
    M77_heading_iNeuron{iNeuron}.tvec(NaNindex) = NaN;
end
for iNeuron = 1:length(M77_neuronsToUse)    % plot the Heading data with overlay of the restricted part for each neuron's PFD
    clf
    plot(M77_heading_degreesTSD.tvec, M77_heading_degreesTSD.data, '.'); hold on
    %     plot(M77_heading_degreesTSD.tvec(M77_headings_to_use{iNeuron}), M77_heading_degreesTSD.data(M77_headings_to_use{iNeuron}), '.');
    plot(M77_heading_iNeuron{iNeuron}.tvec, M77_heading_iNeuron{iNeuron}.data, '.')
    title(num2str(iNeuron))   % how much of the TC is included = 1 - fractiontouse. If fractiontouse = 0.2, then data is restricted to that 80% of TC.
    % If full tuning curve with is 180 deg, then a 0.2 value would mean that we are keeping .8(360-180/360)  40% of the overall session data.
    xlabel('time')
    ylabel('Heading')
    pause
end

