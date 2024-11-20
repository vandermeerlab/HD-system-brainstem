function [slopes, rvalues, neuronList, numSess, cfg_out] = AHV_pearson_correlation(cfg_in, tfilelist)
%2024-11-18. JJS. Calculate the Pearson correlation for each side of the AHV tuning curve (0-90deg/s), as outlined in Graham et al., 2023. This correlation is calculated on the binned 
% (averaged) tuning curve, (not on the raw data).
%  * Might make sense in future to calcaulte and store the binned tuning curves in each folder, so I don't have to recalculate each time. 
% Remember: CW = negative AHV values. CCW = positive AHV values. The AHV xaxis will be CW, CCW from left to right. 

% Inputs: 
%       tfilelist       - a cell array with the dir of each NPH neuron. Like {C:\Jeff\U01\datatouse\M039\M039-2020-08-21-1\M039-2020-08-21-1-TT01_1.t}
%
% Outputs:
%       slopes          - the CW and CCW slopes for the best fit lines of the tuning curve (0-90deg/s). For instance,  0.025 spikes/°/s. 
%                               .CW = clockwise slope. .CCW = counterclockwise slope. .CWcrit = 0 or 1 for meeting minimum slope. .CCWcrit = 0 or 1 for meetinging min slope.
%       rvalues         - Pearson r value for FR vs. AHV for the range from 0 to 90deg/s. 
%                               .CW = r val. for clockwise. .CCW = r val. for counterclockwise. .CWcrit = 0 or 1 for minimum r value. .CCWcrit = 0 or 1 for min value. 
%       neuronList      - easier to read list of neurons that were used in this analysis

cfg_def.min_slope = 0.025;   % from Graham et al.
cfg_def.min_r_value = 0.5;   % from Graham et al.
cfg_out = ProcessConfig2(cfg_def, cfg_in);
slopes = [];
rvalues = []; 
doPlot = 1;
doSave = 1; 

if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
for iNeuron = 1:length(tfilelist)
    [path, neuron_to_use, ext] = fileparts(tfilelist{iNeuron});
    neuronList{iNeuron} = neuron_to_use; 
    disp(neuron_to_use)
    if strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder
        pushdir(path);
        SSN = HD_GetSSN;
        EvalKeys
        numSess = numSess + 1;
    end
    neuronID = strcat(neuron_to_use, ext);
    cfg_spikes.fc = {neuronID};  % do one neuron at a time
    myCell = LoadSpikesJeff(cfg_spikes);
    cellname{iNeuron} = neuronID; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.)
    
    %% get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.minOcc = 100; % the units for this are in samples. Sampling rate for rotation encoder = 200Hz. 1 sample = 5ms. This is the cutoff in Jalina's paper. They accepted anything with at least 0.5 sec data.
    cfg_AHV.subsample_factor = 10;
    cfg_AHV.nBins = 69; % 6 degree/s bins from -204 to 204 works out to 6 degree/s bins. Same as Taube lab papers. From -90 to +90 is 31 bins. 
    cfg_AHV.binEdges = {linspace(-204, 204, cfg_AHV.nBins)};  % cfg_AHV.binEdges(36) == 0. 
    
    [AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, myCell); % position of 'cfg' and 'S' got reversed at some point
    
    x = tc_out.binCenters; x = x';  % tc_out.binCenters(34) is the innermost negative bin [-6 to 0]. tc_out.binCenters(35) is the innermost positive bin [0 to +6].
    y = tc_out.tc;    y = y';
    idx = isnan(y);
    x(:,2) = ones(length(x),1);
    
    % Corr for Positive AHV values
    posIndex = tc_out.binCenters > 0 & tc_out.binCenters <= 90;
    posBins = tc_out.binCenters(posIndex);
    disp(num2str(posBins))
    [bPos,~,~,~,statsPos] = regress(y(posIndex),x(posIndex,:));  % b is the slope (spk/s/deg/s). stats is Rsq, F, p-value, & error variance. Rsq is the coefficient of determination.
    X.bPos(iNeuron) = bPos(1);
    X.rsqPos(iNeuron) = statsPos(1);
    X.pPos(iNeuron) = statsPos(3);
    temp2 = corrcoef(x(posIndex' & ~idx,1),y(posIndex' & ~idx));
    X.Rpos(iNeuron) = temp2(1,2);   % this is the r value (correlation coefficient)
    
    %% Corr for Negative AHV values
    negIndex = tc_out.binCenters < 0 & tc_out.binCenters >=-90;
    negBins = tc_out.binCenters(negIndex);
    disp(num2str(negBins))
    [bNeg,~,~,~,statsNeg] = regress(y(negIndex),x(negIndex,:));
    X.bNeg(iNeuron) = bNeg(1);
    X.rsqNeg(iNeuron) = statsNeg(1);
    X.pNeg(iNeuron) = statsNeg(3);
    temp3 = corrcoef(x(negIndex' & ~idx,1),y(negIndex' & ~idx));
    X.Rneg(iNeuron) = temp3(1,2);
    
    if doPlot 
    plot(tc_out.binCenters, tc_out.tc(iNeuron,:), '.', 'MarkerSize', 25); hold on
    
    
    
    
    
    
    
    
    
    end
    if doSave
        
        
    end
    
    
    
    
    
    
    
end