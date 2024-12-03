function [ipsi_shuff, contra_shuff, neuronList, IC, cfg_out] = calc_circ_shift_AHV_IPSI_CONTRA(cfg_in, tfilelist, X)
%2024-11-22. JJS. This function implements the shuffle procedure for AHV tuning classification, roughly as described in Graham et al., 2023.
%   In Jalina's paper, she adds about 5 sec. to each spike train (relative to the angular position data). Here, I am doing a circular shift with randomly selected times.
% Inputs:
%           tfilelist - cell array list of NPH neurons for analysis. Save in same dir as session data.
%           X         - structure with fields that contain info about the slopes and r values of the true AHV tuning curves. Calculated from AHV_pearson_correlation.m
%                       X.ipsi.r, X.contra.r  - correlation values (r) for the ispi and contra sides (-90 to +90 deg/s) of the tuning curves
%                                          Output X should be obtained with this function   AHV_pearson_correlation_CONTRA_IPSI.m
% Outputs:
%           ipsi_shuff	- shuffle distribution of r values for circ.-shifted spike train (ipsi half of TC)
%           IC.r_sign_ipsi 	- direction of slope. -1 = negative slope. +1 = positive slope.
%           IC.rank_ipsi    - what perentile rank each neuron is in terms of the distribution of shuffled r values. Values closer to 100% indicate a better fit and unlikely due to chance.
%           IC.ipsi_pass   	- whether the neuron passes the r-value shuffle test. +1 = neuron has passed (>95th percentile)
%           IC.per_ipsi_pass- percentage of neurons where the ipsi half of the AHV TC passes threshold (>95th percentile)
%           IC.any_pass     - whether the neuron passes r value shuffle on EITHER the ipsi side or the contra side 
%           IC.per_any_pass - percentage of neurons  where EITHER the ipsi side or the contra side passes threshold
%           ipsi_hem        - this tells you which turn direction is considered ipsilateral (for this neuron, given which hemisphere it was recorded from)

cfg_def.doPlot = 1;
cfg_def.rank_criterion = .05;  % percentage cutoff. 0.025 = 2.5% (bottom 2.5% and top 2.5%)
cfg_def.numShuff = 1000;
cfg_def.startNeuron = 1;
cfg_def.endNeuron = length(tfilelist);
cfg_out = ProcessConfig2(cfg_def, cfg_in);
cfg_out.functionUsed = 'calc_circ_shift_AHV_IPSI_CONTRA.m';

pass_index = round(cfg_out.rank_criterion * cfg_out.numShuff);    % number of samples of the shuffle that the true value needs to be greater(or less) than to pass criterion

if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
for iNeuron = cfg_out.startNeuron : cfg_out.endNeuron
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
    
    %% get True AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.minOcc = 100; % the units for this are in samples. Sampling rate for rotation encoder = 200Hz. 1 sample = 5ms. This is the cutoff in Jalina's paper. They accepted anything with at least 0.5 sec data.
    cfg_AHV.subsample_factor = 10;
    cfg_AHV.nBins = 69; % 6 degree/s bins from -204 to 204 works out to 6 degree/s bins. Same as Taube lab papers. From -90 to +90 is 31 bins.
    cfg_AHV.binEdges = {linspace(-204, 204, cfg_AHV.nBins)};  % cfg_AHV.binEdges(36) == 0.
    cfg_AHV.binSize = median(diff(cfg_AHV.binEdges{1}));
    
%     [~, tc_out] = AHV_tuning(cfg_AHV, myCell); % position of 'cfg' and 'S' got reversed at some point
    
    [endtimetouse, csc_tsd] = findSessionEndTime;
    tvec = csc_tsd.tvec;
    %% Determine Hemisphere
    EvalKeys;
    fc = FindFiles('*.t');
    [a, b, c] = fileparts(fc);
    d = strcat(b,c);
    temp = strcmp(neuronID, d);
    assert(sum(temp) == 1)
    index = find(temp);
    HemisphereToUse = ExpKeys.Hemisphere{index};
    
    %% Caclulate the Circularly Shifted Spike Train
    for iShuff = 1: cfg_out.numShuff
        h = waitbar(iShuff/cfg_out.numShuff);
        disp(num2str(iShuff));
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
        
        %% Correlation for Positive AHV values (CCW)
        CCWindex = tc_out.binCenters > 0 & tc_out.binCenters <= 90;
        tempCCW = corrcoef(x(CCWindex' & ~idx,1),y(CCWindex' & ~idx));
        CCW_shuff{iNeuron}(iShuff) = tempCCW(1,2);   % this is the r value (correlation coefficient)
        
        %% Correlation for Negative AHV values  (CW)
        CWindex = tc_out.binCenters < 0 & tc_out.binCenters >=-90;
        tempCW = corrcoef(x(CWindex' & ~idx,1),y(CWindex' & ~idx));
        CW_shuff{iNeuron}(iShuff) = tempCW(1,2);
        
        % Translate to IPSI/CONTRA
        if strcmp(HemisphereToUse, 'L')   % Leftwards = CCW. Rightwards = CW. L hem. neuron = CCW. 
            ipsi_shuff{iNeuron} = CCW_shuff{iNeuron};
            contra_shuff{iNeuron} = CW_shuff{iNeuron};
            ipsi_hem = 'CCW'; 
        elseif strcmp(HemisphereToUse, 'R')  % R hem. neuron = CW. 
            ipsi_shuff{iNeuron} = CW_shuff{iNeuron};
            contra_shuff{iNeuron} = CCW_shuff{iNeuron};
            ipsi_hem = 'CW'; 
        else
            warning('hemisphere info not correct')
            ipsi_shuff{iNeuron} = {};
            contra_shuff{iNeuron} = {};
        end
        r_ipsi_shuff_sorted{iNeuron} = sort(ipsi_shuff{iNeuron});
        r_contra_shuff_sorted{iNeuron} = sort(contra_shuff{iNeuron});
    end
end
total_neurons = length(tfilelist);
%% Preallocate
IC.r_ipsi_true = X.ipsi.r;
IC.r_contra_true = X.contra.r;

IC.r_ipsi_shuff = ipsi_shuff;
IC.r_contra_shuff = contra_shuff;

IC.r_ipsi_shuff_sorted = r_ipsi_shuff_sorted;
IC.r_contra_shuff_sorted = r_contra_shuff_sorted;

IC.ipsi_pass = zeros(1,cfg_out.endNeuron);
IC.contra_pass = zeros(1,cfg_out.endNeuron);

IC.any_pass = nan(1,cfg_out.endNeuron);

IC.r_sign_ipsi = nan(1,cfg_out.endNeuron);
IC.r_sign_contra = nan(1,cfg_out.endNeuron);

ix_ipsi_neg = nan(1,cfg_out.endNeuron);
ix_ipsi_pos = nan(1,cfg_out.endNeuron);

ix_contra_neg = nan(1,cfg_out.endNeuron);
ix_contra_pos = nan(1,cfg_out.endNeuron);

IC.low_index = pass_index;
IC.high_index = cfg_out.numShuff - pass_index;

%% Determine if r Criterion met and Rank
for iNeuron = cfg_out.startNeuron : cfg_out.endNeuron
    % IPSI loop  
    % negative slope. true r value must be lower than criterion index (place/rank)
    if X.ipsi.r(iNeuron) < 0   % r value negative (inverse relationship btwn FR and AHV)
        IC.r_sign_ipsi(iNeuron) = -1;
        crit = r_ipsi_shuff_sorted{iNeuron}(IC.low_index);
        [~, ix_ipsi_neg(iNeuron)] = min(abs(r_ipsi_shuff_sorted{iNeuron} - X.ipsi.r(iNeuron)));
        if X.ipsi.r(iNeuron) <= crit  % threshold is rank lower than crit. For instance, 50th place or earlier (negative side of r-value distribution)
            IC.ipsi_pass(iNeuron) = 1;
        end
        % positive slope. true r value must be greater than criterion index (place/rank)
    elseif X.ipsi.r(iNeuron) > 0  % r value positive (direct relationship btwn FR and AHV)
        IC.r_sign_ipsi(iNeuron) = 1;
        crit = r_ipsi_shuff_sorted{iNeuron}(IC.high_index);
        [~, ix_ipsi_pos(iNeuron)] = min(abs(r_ipsi_shuff_sorted{iNeuron} - X.ipsi.r(iNeuron)));
        if X.ipsi.r(iNeuron) >= crit % threshold is rank higher than crit. For instance, 950th place or later (positive side of r-value distribution)
            IC.ipsi_pass(iNeuron) = 1;
        end
    else
        error('slope ipsi exactly equals zero')
    end
    % CONTRA loop 
    % negative slope. true r value must be lower than criterion index (place/rank)
    if X.contra.r(iNeuron) < 0
        IC.r_sign_contra(iNeuron) = -1;
        crit = r_contra_shuff_sorted{iNeuron}(IC.low_index);
        [~, ix_contra_neg(iNeuron)] = min(abs(r_contra_shuff_sorted{iNeuron} - X.contra.r(iNeuron)));
        if X.contra.r(iNeuron) <= crit
            IC.contra_pass(iNeuron) = 1;
        end
        % positive slope. true r value must be greater than criterion index (place/rank)
    elseif X.contra.r(iNeuron) > 0
        IC.r_sign_contra(iNeuron) = 1;
        crit = r_contra_shuff_sorted{iNeuron}(IC.high_index);
        [~, ix_contra_pos(iNeuron)] = min(abs(r_contra_shuff_sorted{iNeuron} - X.contra.r(iNeuron)));
        if X.contra.r(iNeuron) >= crit
            IC.contra_pass(iNeuron) = 1;
        end
    else
        error('slope contra exactly equals zero')
    end
    if IC.ipsi_pass(iNeuron) || IC.contra_pass(iNeuron)
        IC.any_pass(iNeuron) = 1;
    end
end
%% Convert to 500-1000 scale
ix_ipsi_neg = repmat(cfg_out.numShuff, 1, total_neurons) - ix_ipsi_neg;   % CW. For negative r values, the true rvalue needs to be in the bottom 5% of shuffle. To convert to a standard percentile, we subtract 1000 - n. If index = 50, then percentile = 95% (not 5%)
ix_contra_neg = repmat(cfg_out.numShuff, 1, total_neurons) - ix_contra_neg; % CCW.

% IPSI
ix_ipsi(IC.r_sign_ipsi == -1) = ix_ipsi_neg(IC.r_sign_ipsi == -1);   % add in the negative r value TC half rankings
ix_ipsi(IC.r_sign_ipsi == 1)  = ix_ipsi_pos(IC.r_sign_ipsi == 1);   % add in the positive r value TC half rankings
IC.rank_ipsi = ix_ipsi./10;

% CONTRA
ix_contra(IC.r_sign_contra == -1) = ix_contra_neg(IC.r_sign_contra == -1);
ix_contra(IC.r_sign_contra == 1)  = ix_contra_pos(IC.r_sign_contra == 1);
IC.rank_contra = ix_contra./10;

IC.per_ipsi_pass = sum(IC.ipsi_pass)/total_neurons;
IC.per_contra_pass = sum(IC.contra_pass)/total_neurons;
IC.per_any_pass = nansum(IC.any_pass)/total_neurons;

IC.ipsi_hem = ipsi_hem;
end

