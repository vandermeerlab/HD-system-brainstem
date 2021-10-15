function [FRxBinT, FRxBinN, FRxBinTsmooth, FRxBinNsmooth, FRxBinTnorm, FRxBinNnorm, TnormSmooth, NnormSmooth, outputIT, binCenters, cfg, cellID, cellname] = makeSaccadeHeatPlot(cfg_in, varargin)
% [] = makeSaccadeHeatPlot
% JJS. 4/2021.

% doPlot = 0;
FRwindow = [-.2 .2];
dt = 0.01;
pethBins = linspace(FRwindow(1), FRwindow(2), diff(FRwindow)/dt+1);
cfg_def.threshH = 10;
cfg_def.threshL = -10;
cfg = ProcessConfig(cfg_def,cfg_in);

process_varargin(varargin);
cellCounterToUse = 0;
cellCounterIndex = 0;

fd = FindFiles('*keys.m');
% FRxBinT = nan(length(fd), length(pethBins));
% FRxBinN = nan(length(fd), length(pethBins));
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN;
    disp(SSN)
    [S] = LoadSpikesJeff;
    if exist(strcat(SSN, '-VT1_proc.mat'))
        %         [S] = LoadSpikesJeff;
        [temporalSaccades, nasalSaccades, ~, ~, ~, ~, ~, ~, ~] = processPupilData2(cfg_in, 'doPlotEverything', 0, 'doPlotThresholds', 0);
        spikefiles = FindFiles('*.t');
        for iCell = 1:length(S.t)
            [a, b, c] = fileparts(spikefiles{iCell});
            cellID = b; 
            cellCounterToUse = cellCounterToUse + 1;   % this counts eye tracking neuron. Neurons from sessions that are skipped (above) are not recorded. This total will = total neurons with eye tracking.
            cellCounterIndex = cellCounterIndex +1;    % this counts every neuron, whether or not it has eye tracking data associated with it. 
            cellID(cellCounterToUse) = cellCounterIndex; % this is the order of the eye tracking neurons, with respect to the total group of neurons. 
            cellname{cellCounterToUse} = S.label{iCell}; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.) 
            myCell = SelectTS([], S, iCell);
            cfg_in = [];
            cfg_in.window = FRwindow;
            cfg_in.dt = dt;
            figure
            [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm(cfg_in, myCell, temporalSaccades);
            subplot(2,1,1); title(cellID)
            subplot(2,1,2); title('Temporal Saccades')
            assert(isequal(pethBins, outputIT))
            mT = histcounts(outputS, outputIT);
            FRxBinT(cellCounterToUse,:) = mT/cfg.dt/length(temporalSaccades);
            clear outputS;
            figure
            [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm(cfg_in, myCell, nasalSaccades);
            subplot(2,1,1); title(cellID)
            subplot(2,1,2); title('Nasal Saccades')
            mN = histcounts(outputS, outputIT);
            FRxBinN(cellCounterToUse,:) = mN/cfg.dt/length(nasalSaccades);
        end
    else
        for iCell = 1:length(S.t)
            cellCounterIndex = cellCounterIndex +1;
%             cellCounterToUse = cellCounterToUse + 1;
%             FRxBinT(cellCounterToUse,:) = nan(1,length(pethBins)-1);
%             FRxBinN(cellCounterToUse,:) = nan(1,length(pethBins)-1);
        end
        disp('no eye tracking file detected. Skipping Session.')
    end
    popdir;
end
A = outputIT;
binCenters = (A(:,1:end-1) + A(:,2:end)) / 2;
%% Temporal Saccades
[FRxBinTsmooth, smoothwindow_a] = smoothdata(FRxBinT,2);
[mT, iT] = max(FRxBinT,[], 2);
FRxBinTnorm = FRxBinT./mT;       % normalized by FR
[bT, sT] = sort(iT);             % sT = sorted, unnormalized
[TnormSmooth, smoothwindow_b] = smoothdata(FRxBinTnorm,2);

%% Nasal Saccades
[FRxBinNsmooth, smoothwindow_c] = smoothdata(FRxBinN,2);
[mN, iN] = max(FRxBinN,[], 2);
FRxBinNnorm = FRxBinN./mN;       % normalized by FR
[bN, sN] = sort(iN);             % sT = sorted, unnormalized
[NnormSmooth, smoothwindow_d] = smoothdata(FRxBinNnorm,2);

% %% Arrange the data in pairs so that you can see temporal and nasal saccades for the same cell, one above the other
% pairs = nan(size(FRxBinT,1), size(FRxBinT,2));
% counter = 0;
% for iPosition = 2:2:2*(size(FRxBinN,1))   % nasal
%     counter = counter + 1;
%     pairs(iPosition, :) = FRxBinN(counter,:);
% end
% counter = 0;
% for iPosition = 1:2:2*(size(FRxBinT,1))-1    % temporal
%     counter = counter + 1;
%     pairs(iPosition, :) = FRxBinT(counter,:);
% end
% 
% pairsNorm = nan(size(FRxBinT,1), size(FRxBinT,2));
% counter = 0;
% for iPosition = 2:2:2*(size(FRxBinN,1))   % nasal
%     counter = counter + 1;
%     pairsNorm(iPosition, :) = NnormSmooth(counter,:);
% end
% counter = 0;
% for iPosition = 1:2:2*(size(FRxBinT,1))-1    % temporal
%     counter = counter + 1;
%     pairsNorm(iPosition, :) = TnormSmooth(counter,:);
% end
% 
% 
% if doPlot == 1
%     
% end





