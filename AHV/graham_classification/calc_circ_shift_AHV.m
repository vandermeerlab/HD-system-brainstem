function [Z, cfg_out] = calc_circ_shift_AHV(cfg_in, tfilelist, X)
%2024-11-22. JJS. This function implements the shuffle procedure for AHV tuning classification, roughly as described in Graham et al., 2023.
%   In Jalina's paper, she adds about 5 sec. to each spike train (relative to the angular position data). Here, I am doing a circular shift with randomly selected times.
% Inputs:
%           tfilelist - cell array list of NPH neurons for analysis. Save in same dir as session data.
%           X         - structure with fields that contain info about the slopes and r values of the true AHV tuning curves
%                       X.Rpos. X.Rneg.  - correlation values (r) for the CCW and CW sides (-90 to +90 deg/s) of the tuning curves
%                                          Output X should be obtained with this function   AHV_pearson_correlation.m
% Outputs:
%           Z         - structure with fields for the rank of each neuron in the shuffle distribution and the percentage of cells that passed the r - criterion
tic
clear Z
cfg_def.rank_criterion = .05;  % percentage. what percent of the distribution the true value needs to exceed to meet threshold.
% For neg. slopes, look at bottom 5%. For pos. slopes, look at top 5%.
cfg_def.numShuff = 1000;
cfg_def.startNeuron = 1;
cfg_def.endNeuron = length(tfilelist);
cfg_out = ProcessConfig2(cfg_def, cfg_in);
total_neurons = cfg_out.endNeuron - cfg_out.startNeuron + 1;

pass_index = round(cfg_out.rank_criterion * cfg_out.numShuff);    % number of samples of the shuffle that the true value needs to be greater(or less) than to pass criterion

r_CCW_true = X.Rpos;  % change the naming to CW & CCW, which it should have been in the first place. 'pos' here refers to positive AHV values. but having AHV, slopes, and r values all be pos/neg is confusing
r_CW_true = X.Rneg;

if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
for iNeuron = cfg_out.startNeuron : cfg_out.endNeuron
    tic
    [path, neuron_to_use, ext] = fileparts(tfilelist{iNeuron});
    neuronList{iNeuron} = neuron_to_use;
    disp(neuron_to_use)
    if strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder
        pushdir(path);
        SSN = HD_GetSSN;
        EvalKeys
    end
    neuronID = strcat(neuron_to_use, ext);
    cfg_spikes.fc = {neuronID};  % do one neuron at a time
    myCell = LoadSpikesJeff(cfg_spikes);
    cellname{iNeuron} = neuronID; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.)
    
%     %% get the true AHV Tuning Curve   % this is already supplied in the variable X (pre-calculated) 
%     cfg_AHV = [];
%     cfg_AHV.minOcc = 100; % the units for this are in samples. Sampling rate for rotation encoder = 200Hz. 1 sample = 5ms. This is the cutoff in Jalina's paper. They accepted anything with at least 0.5 sec data.
%     cfg_AHV.subsample_factor = 10;
%     cfg_AHV.nBins = 69; % 6 degree/s bins from -204 to 204 works out to 6 degree/s bins. Same as Taube lab papers. From -90 to +90 is 31 bins.
%     cfg_AHV.binEdges = {linspace(-204, 204, cfg_AHV.nBins)};  % cfg_AHV.binEdges(36) == 0.
%     binSize = median(diff(cfg_AHV.binEdges{1}));
%     
%     [~, tc_out] = AHV_tuning(cfg_AHV, myCell); % position of 'cfg' and 'S' got reversed at some point
    
    [endtimetouse, csc_tsd] = findSessionEndTime;
    tvec = csc_tsd.tvec;
    
    %% Caclulate the Circularly Shifted Spike Train
    for iShuff = 1: cfg_out.numShuff
        disp(num2str(iShuff))
        h = waitbar(iShuff/cfg_out.numShuff);
        %     disp(iShuff);
        r(iShuff) = randsample(tvec,1); % choose a random value from the session times tvec
        
        spike_shift = myCell.t{1} + r(iShuff);
        
        % SUBRACT SESSION DURATION FROM VALUES GREATER THAN THE LAST SESSION TIMESTAMP
        greater_list = spike_shift > endtimetouse;
        
        if ~isempty(greater_list)
            spikes_to_stay = spike_shift <= endtimetouse;
            t1 = spike_shift(spikes_to_stay);
            t2 = spike_shift(greater_list) - endtimetouse;
            spike_shift_new = vertcat(t2, t1);
            num_new = length(spike_shift_new);
            num_old = length(myCell.t{1});
            assert(num_new == num_old);
        else
            disp('warning, spike train was the same after shift')
            spike_shift_new = spike_shift;
        end
        myCellShuff = ts({spike_shift_new});
        
        [~, tc_out] = AHV_tuning(cfg_AHV, myCellShuff); % position of 'cfg' and 'S' got reversed at some point
        x = tc_out.binCenters; x = x';  % tc_out.binCenters(34) is the innermost negative bin [-6 to 0]. tc_out.binCenters(35) is the innermost positive bin [0 to +6].
        y = tc_out.tc;    y = y';
        idx = isnan(y);
        x(:,2) = ones(length(x),1);
        
        % Corr for Positive AHV values (CCW)
        CCW_index = tc_out.binCenters > 0 & tc_out.binCenters <= 90;
        tempPos = corrcoef(x(CCW_index' & ~idx,1),y(CCW_index' & ~idx));
        r_CCW_shuff{iNeuron}(iShuff) = tempPos(1,2);   % this is the r value (correlation coefficient)
        
        %% Corr for Negative AHV values (CW)
        CW_index = tc_out.binCenters < 0 & tc_out.binCenters >=-90;
        tempNeg = corrcoef(x(CW_index' & ~idx,1),y(CW_index' & ~idx));
        r_CW_shuff{iNeuron}(iShuff) = tempNeg(1,2);
    end
    r_CCW_shuff_sorted{iNeuron} = sort(r_CCW_shuff{iNeuron});
    r_CW_shuff_sorted{iNeuron} = sort(r_CW_shuff{iNeuron});
    toc
end
toc

%% Determine whether neurons meet the threshold
Z.r_CCW_true = r_CCW_true; 
Z.r_CW_true = r_CW_true; 
Z.r_CCW_shuff = r_CCW_shuff;
Z.r_CW_shuff = r_CW_shuff;
Z.r_CCW_shuff_sorted = r_CCW_shuff_sorted; 
Z.r_CW_shuff_sorted = r_CW_shuff_sorted;

Z.CW_pass = zeros(1,cfg_out.endNeuron);
Z.CCW_pass = zeros(1,cfg_out.endNeuron);
Z.any_pass = nan(1,cfg_out.endNeuron);
Z.r_sign_CW = nan(1,cfg_out.endNeuron);
Z.r_sign_CCW = nan(1,cfg_out.endNeuron);
ix_CW = nan(1,cfg_out.endNeuron);
ix_CW_neg = nan(1,cfg_out.endNeuron);
ix_CW_pos = nan(1,cfg_out.endNeuron);
ix_CCW_neg = nan(1,cfg_out.endNeuron);
ix_CCW_pos = nan(1,cfg_out.endNeuron);

Z.low_index = pass_index;
Z.high_index = cfg_out.numShuff - pass_index;

%% Determine if r Criterion met and Rank
for iNeuron = cfg_out.startNeuron : cfg_out.endNeuron
    % Negative AHV loop  (CW)
    if r_CW_true(iNeuron) < 0   % r value negative (inverse relationship btwn FR and AHV)
        Z.r_sign_CW(iNeuron) = -1;
        crit = r_CW_shuff_sorted{iNeuron}(Z.low_index);
        [~, ix_CW_neg(iNeuron)] = min(abs(r_CW_shuff_sorted{iNeuron} - r_CW_true(iNeuron)));
        if r_CW_true(iNeuron) <= crit
            Z.CW_pass(iNeuron) = 1;
        end
    elseif r_CW_true(iNeuron) > 0  % r value positive (direct relationship btwn FR and AHV)
        Z.r_sign_CW(iNeuron) = 1;
        crit = r_CW_shuff_sorted{iNeuron}(Z.high_index);
        [~, ix_CW_pos(iNeuron)] = min(abs(r_CW_shuff_sorted{iNeuron} - r_CW_true(iNeuron)));
        if r_CW_true(iNeuron) >= crit
            Z.CW_pass(iNeuron) = 1;
        end
    else
        error('slope_CW exactly equals zero')
    end
    % Positive AHV loop  (CCW)
    if r_CCW_true(iNeuron) < 0
        Z.r_sign_CCW(iNeuron) = -1;
        crit = r_CCW_shuff_sorted{iNeuron}(Z.low_index);
        [~, ix_CCW_neg(iNeuron)] = min(abs(r_CW_shuff_sorted{iNeuron} - r_CCW_true(iNeuron)));
        if r_CCW_true(iNeuron) <= crit
            Z.CCW_pass(iNeuron) = 1;
        end
    elseif r_CCW_true(iNeuron) > 0
        Z.r_sign_CCW(iNeuron) = 1;
        crit = r_CCW_shuff_sorted{iNeuron}(Z.high_index);
        [~, ix_CCW_pos(iNeuron)] = min(abs(r_CCW_shuff_sorted{iNeuron} - r_CCW_true(iNeuron)));
        if r_CCW_true(iNeuron) >= crit
            Z.CCW_pass(iNeuron) = 1;
        end
    else
        error('slope_CCW exactly equals zero')
    end
    if Z.CW_pass(iNeuron) || Z.CCW_pass(iNeuron)
        Z.any_pass(iNeuron) = 1;
    end
end
%% Convert to 500-1000 scale
ix_CWneg = repmat(cfg_out.numShuff, 1, total_neurons) - ix_CW_neg;   % CW. For negative r values, the true rvalue needs to be in the bottom 5% of shuffle. To convert to a standard percentile, we subtract 1000 - n. If index = 50, then percentile = 95% (not 5%)
ix_CCWneg = repmat(cfg_out.numShuff, 1, total_neurons) - ix_CCW_neg; % CCW. 

% CW
ix_CW(Z.r_sign_CW == -1) = ix_CWneg(Z.r_sign_CW == -1);   % add in the negative r value TC half rankings 
ix_CW(Z.r_sign_CW == 1)  = ix_CW_pos(Z.r_sign_CW == 1);   % add in the positive r value TC half rankings 
Z.rank_CW = ix_CW./10;

% CCW
ix_CCW(Z.r_sign_CCW == -1) = ix_CCWneg(Z.r_sign_CCW == -1);
ix_CCW(Z.r_sign_CCW == 1)  = ix_CCW_pos(Z.r_sign_CCW == 1);
Z.rank_CCW = ix_CCW./10;

Z.per_neg_pass = sum(Z.CW_pass)/total_neurons;
Z.per_pos_pass = sum(Z.CCW_pass)/total_neurons;
Z.per_any_pass = nansum(Z.any_pass)/total_neurons;
end








