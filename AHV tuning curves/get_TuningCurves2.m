function [TC_all, TC_norm, TC_cellID] = get_TuningCurves2(cfg_in, tfilelist)
% JJS. 2024-10-01. This function calculated the AHV tuning curve for all neurons in the tfilelist and plots them, sorted, as an imagesc image.

% Inputs:

% Outputs:
%           sessCounter: how many sessions the function went in to to generate the data. This should equal 109 for NPH neurons (as of 2024-10-01).

subtractStartTime = 1;
cfg_def = [];
cfg = ProcessConfig(cfg_def,cfg_in);

sessCounter = 0;
TC_all.tc = [];
TC_all.bins = [];
TC_cellID = [];

if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
for iNeuron = 1:length(tfilelist)
    [path, neuron_to_use, ~] = fileparts(tfilelist{iNeuron});
    if strcmp(pwd, path) == 0
        pushdir(path);
        SSN = HD_GetSSN; disp(SSN)
        sessCounter = sessCounter + 1;
        EvalKeys
        tS = FindFiles('*.t'); % list of neurons w/ path
        [~, neuronIDs, ~] = fileparts(tS);
        %% Load Spikes
        cfg_S = [];
        cfg_S.uint = '64';
        S = LoadSpikes(cfg_S);
        if subtractStartTime == 1 % New cheetah versions have timestamps that start at UNIX time zero
            events_ts = LoadEvents([]);
            assert(strcmp(events_ts.label{1}, 'Starting Recording')) = 1;
            for iC = 1:length(S.t)
                S.t{iC} = S.t{iC} - events_ts.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
            end
        end
        %% Calculate the AHV tuning curve
        cfg_AHV = []; % default number of bins here is 100
        cfg_AHV.subsample_factor = 10;
        [~, tc_out] = AHV_tuning(cfg_AHV, S);
    end
    index = strcmp({neuron_to_use}, neuronIDs);  % this result should be a nCell x 1 logical. if [0 1 0] for instance, then the second neuron in S is our spike train to use.
    rows_to_use = find(index);
    TCs_to_use = tc_out.tc(rows_to_use,:);
    
    %% Figure out which tuning curves are NPH and therefore which to concatenate
    if ~isempty(rows_to_use)
        numCells = size(tc_out.tc, 1);
        TC_all.tc = vertcat(TC_all.tc, TCs_to_use);
        TC_all.bins = vertcat(TC_all.bins, repmat(tc_out.binCenters, rows_to_use, 1));
        for iLabel = 1:length(rows_to_use)
            TC_cellID = vertcat(TC_cellID, {S.label{rows_to_use(iLabel)}});
        end
    end
    popdir;
end

TC_norm.tc = [];
TC_norm.bins = TC_all.bins;
for iCell = 1:size(TC_all.tc,1)
    TC_norm.tc(iCell,:) = TC_all.tc(iCell,:)./max(TC_all.tc(iCell,:));
end

[M, I] = max(TC_all.tc');
[B, isort] = sort(I);
imagesc(TC_all.bins(1,:), 1:length(TC_norm.tc), TC_norm.tc(isort, :))
cmap = parula(256);
% Make lowest one white
cmap(1,:) = 0;
colormap(cmap);


% %% plot it 
% [M, I] = max(TC_all.tc');
% 
% [B, isort] = sort(I);
% 
% pcolor(TC_all.bins(1,:), 1:length(TC_norm.tc), TC_norm.tc(isort, :))
% 
% cmap = jet(256);
% % Make lowest one black
% cmap(1,:) = 0;
% colormap(cmap);


