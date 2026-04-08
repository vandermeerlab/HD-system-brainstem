function [Z] = makeSaccadeHeatPlotAllneurons_temporal_nasal(cfg_in, tfilelist)
% [] = makeSaccadeHeatPlot
% JJS. 4/2021.
FontSize = 25;
LineWidth = 3;
doPlot = 0;
cfg_def.doPlot = 1;
cfg_def.FRwindow = [-.2 .2];
cfg_def.dt = 0.005;
cfg_def.threshH = 10;
cfg_def.threshL = -10;
cfg = ProcessConfig(cfg_def,cfg_in);
pethBins = linspace(cfg.FRwindow(1), cfg.FRwindow(2), diff(cfg.FRwindow)/cfg.dt+1);
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
    
    %     %% Isolate the spike train from our neuron of interest
    %     tS = FindFiles('*.t');
    %     [a b c] = fileparts(tS);
    %     which_neuron = strcmp(neuron_to_use, b);
    %     neuronIndex = find(which_neuron);
    %     myCell = SelectTS([], S, neuronIndex);
    
    %% Calculate Temporal Saccades
    cfg_in = [];
    cfg_in.window = cfg.FRwindow;
    cfg_in.dt = cfg.dt;
    [outputS, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, temporalSaccades, 'doPlot', doPlot);
    assert(isequal(pethBins, outputIT))
    mT = histcounts(outputS, outputIT);
    FRxBinT(iNeuron,:) = mT/cfg.dt/length(temporalSaccades);
    
    %% Calculate Nasal Saccades
    [outputS, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, nasalSaccades, 'doPlot', doPlot);
    mN = histcounts(outputS, outputIT);
    FRxBinN(iNeuron,:) = mN/cfg.dt/length(nasalSaccades);
    popdir;
end
A = outputIT;
binCenters = (A(:,1:end-1) + A(:,2:end)) / 2;
%% Temporal Saccades
[FRxBinTsmooth, ~] = smoothdata(FRxBinT,2);
[mT, iT] = max(FRxBinT,[], 2);
FRxBinTnorm = FRxBinT./mT;       % normalized by FR
[~, sT] = sort(iT);             % sT = sorted, unnormalized
[TnormSmooth, ~] = smoothdata(FRxBinTnorm,2);

%% Nasal Saccades
[FRxBinNsmooth, ~] = smoothdata(FRxBinN,2);
[mN, iN] = max(FRxBinN,[], 2);
FRxBinNnorm = FRxBinN./mN;       % normalized by FR
[~, sN] = sort(iN);             % sT = sorted, unnormalized
[NnormSmooth, ~] = smoothdata(FRxBinNnorm,2);

Z.FRxBinT = FRxBinT;
Z.FRxBinN = FRxBinN;
Z.FRxBinTsmooth = FRxBinTsmooth;
Z.FRxBinNsmooth = FRxBinNsmooth;
Z.FRxBinTnorm = FRxBinTnorm;
Z.FRxBinNnorm = FRxBinNnorm;
Z.TnormSmooth = TnormSmooth;
Z.NnormSmooth = NnormSmooth;
Z.outputIT = outputIT;
Z.binCenters = binCenters;
Z.cfg = cfg;
Z.cellname = cellname;

Z.sT = sT;
Z.sN = sN; 

%% Plot It 
if cfg.doPlot
    clf
    subplot(1,2,1)
    imagesc(Z.binCenters, 1:length(Z.TnormSmooth), Z.TnormSmooth(Z.sT,:));
    title('Temporal Saccades')
    xlabel('time (ms)')
    ylabel('neuron #')
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'LineStyle', '-', 'Color', 'r', 'LineWidth', LineWidth)
    
    subplot(1,2,2)
    imagesc(Z.binCenters, 1:length(Z.NnormSmooth), Z.NnormSmooth(Z.sN,:));
    title('Nasal Saccades')
    xlabel('time (ms)')
    ylabel('neuron #')
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'LineStyle', '-', 'Color', 'r', 'LineWidth', LineWidth)
end







