function [neuronList, AHV_r, AHV_slope, AHV_index, AHV_score] = calc_AHV_score_tfilelist(cfg_in, tfilelist)
% JJS. 2024-10-17.  This function takes a cell array of neuron paths (tfilelist) and calculates the best line fit (r) and slope (in spk/s/deg/s), and AHV score 
%                   defined as r x slope x 100 for each turn direction (CW & CCW), and taking the larger of those 2 values. AHV values at +-5 deg/s AHV are omitted. 
% Inputs:           cfg_in      - structure with configuration variables
%                   tfilelist   - cell array. each cell is the path for one neuron, like {;'C:\Jeff\U01\datatouse\M039\M039-2020-08-21-1\M039-2020-08-21-1-TT01_1.t'}

% Outputs:          neuronList  - a list of neuron IDs that were analyzed
%                   AHV_r       - 2 x n neuron double of correlation r values. Column 1 is IPSI side. Column 2 is CONTRA side. 
%                   AHV_slope   - 2 x n neuron double of best fit line slopes. Units are spk/s/deg/sec. Column 1 is IPSI side. Column 2 is CONTRA side.
%                   AHV_score   - 2 x n neuron double of correlation AHV scores. As defined by Taube (personal communication). Column 1 is IPSI side. Column 2 is CONTRA side.
%                   AHV_index   - 1 x n neuron double of whichever AHV_score is greater for each neuron.   

%                       rsqAll:     [1×n double]     Rsq statistic for each neuron, fit with a single line
%                       pAll:       [1×n double]     p-values for each neuron, fit with a single line
%                       coeffsAll:  [n×2 double]     1st row = slope (Beta value); 2nd row = y-intercept. Single line fit
%                       rsqPos:     [1×n double]     Rsq for positive part of tuning curve (above zero deg/sec)
%                       pPos:       [1×n double]     p-values for positive part of tuning curve
%                       coeffsPos:  [n×2 double]     coeffs for positive part of tuning curve
%                       rsqNeg:     [1×n double]     Rsq for negative part of tuning curve
%                       pNeg:       [1×n double]     p-values for negative part of tuning curve
%                       coeffsNeg:  [n×1 double]     coeffs for negative part of tuning curve
sessCounter = 0;
if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
for iNeuron = 1:length(tfilelist)
    [path, neuron_to_use, ext] = fileparts(tfilelist{iNeuron});
    disp(neuron_to_use)
    if strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder
        pushdir(path);
        SSN = HD_GetSSN;
        EvalKeys
        sessCounter = sessCounter + 1;
        load(FindFile('*saccades-edited.mat'), 'nasalSaccades', 'temporalSaccades')
        keep = find(~isnan(temporalSaccades)); % get rid of NaNs
        temporalSaccades = temporalSaccades(keep);  %#ok<*FNDSB>
        keep = find(~isnan(nasalSaccades));
        nasalSaccades = nasalSaccades(keep);
    end
    neuronID = strcat(neuron_to_use, ext);
    cfg_spikes.fc = {neuronID};  % do one neuron at a time
    myCell = LoadSpikesJeff(cfg_spikes);
    cellname{iNeuron} = neuronID; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.)

    [S] = LoadSpikesJeff;
    if ~isempty(S.t)
        % get AHV Tuning Curve
        cfg_AHV = [];
        cfg_AHV.subsample_factor = 10;
        [AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, S); % position of 'cfg' and 'S' got reversed at some point
        AHV_dt = median(diff(AHV_tsd.tvec));



















end

