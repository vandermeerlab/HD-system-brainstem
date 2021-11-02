function [FRxBinT, FRxBinN, FRxBinTsmooth, FRxBinNsmooth, FRxBinTnorm, FRxBinNnorm, TnormSmooth, NnormSmooth, outputIT, binCenters, cfg, cellID, cellname] = makeSaccadeHeatPlot(cfg_in, varargin)
% [] = makeSaccadeHeatPlot
% JJS. 4/2021.

doPlot = 0;
cfg_def.FRwindow = [-.2 .2];
cfg_def.dt = 0.005;
cfg_def.threshH = 10;
cfg_def.threshL = -10;
cfg = ProcessConfig(cfg_def,cfg_in);

pethBins = linspace(cfg.FRwindow(1), cfg.FRwindow(2), diff(cfg.FRwindow)/cfg.dt+1);

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
    if strcmp(SSN, 'M293-2021-08-19') ~= 1
        [S] = LoadSpikesJeff;
        if exist(strcat(SSN, '-VT1_proc.mat'))
            %         [S] = LoadSpikesJeff;
            spikefiles = FindFiles('*.t');
            [temporalSaccades, nasalSaccades, combinedSaccades, index_tP_final, index_nP_final, tsdH, tsdV, diffH, diffV] = processPupilData2(cfg);
            for iCell = 1:length(S.t)
                [a, b, c] = fileparts(spikefiles{iCell});
                celltitle = b;
                cellCounterToUse = cellCounterToUse + 1;   % this counts eye tracking neuron. Neurons from sessions that are skipped (above) are not recorded. This total will = total neurons with eye tracking.
                cellCounterIndex = cellCounterIndex +1;    % this counts every neuron, whether or not it has eye tracking data associated with it.
                cellID(cellCounterToUse) = cellCounterIndex; % this is the order of the eye tracking neurons, with respect to the total group of neurons.
                cellname{cellCounterToUse} = S.label{iCell}; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.)
                myCell = SelectTS([], S, iCell);
                cfg_in = [];
                cfg_in.window = cfg.FRwindow;
                cfg_in.dt = cfg.dt;
                if doPlot == 1
                    figure(1)
                    subplot(2,1,1); title(celltitle)
                    subplot(2,1,2); title('Temporal Saccades')
                end
                [outputS, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, temporalSaccades, 'doPlot', doPlot);
                assert(isequal(pethBins, outputIT))
                mT = histcounts(outputS, outputIT);
                FRxBinT(cellCounterToUse,:) = mT/cfg.dt/length(temporalSaccades);
                if doPlot == 1
                    figure(2)
                    subplot(2,1,1); title(celltitle)
                    subplot(2,1,2); title('Nasal Saccades')
                end
                [outputS, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, nasalSaccades, 'doPlot', doPlot);
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





