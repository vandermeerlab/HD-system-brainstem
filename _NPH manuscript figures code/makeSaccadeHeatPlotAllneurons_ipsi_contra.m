function [Z] = makeSaccadeHeatPlotAllneurons_ipsi_contra(cfg_in, tfilelist)
% [] = makeSaccadeHeatPlot
% JJS. 4/2021.
FontSize = 30;
LineWidth = 3;
doPlot = 0;
cfg_def.doPlot = 0;
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
    clear ExpKeys
    if strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder
        pushdir(path);
        SSN = HD_GetSSN;
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
    EvalKeys
    
    %% Isolate the spike train from our neuron of interest
    tS = FindFiles('*.t');
    [a, b, c] = fileparts(tS);
    which_neuron = strcmp(neuron_to_use, b);
    neuronIndex = find(which_neuron);
    
    if strcmp(ExpKeys.Hemisphere{neuronIndex}, 'R')
        ipsi_saccades = nasalSaccades;
        contra_saccades = temporalSaccades;
    elseif strcmp(ExpKeys.Hemisphere{neuronIndex}, 'L')
        contra_saccades = nasalSaccades;
        ipsi_saccades = temporalSaccades;
    else
        error('issue with neuron hemisphere')
    end
    %% Calculate IPSI Saccades
    cfg_in = [];
    cfg_in.doPlot = 0;
    cfg_in.doBar = 0;
    cfg_in.window = cfg.FRwindow;
    cfg_in.dt = cfg.dt;
    [outputIPSI, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, ipsi_saccades, 'doPlot', doPlot);
    assert(isequal(pethBins, outputIT))
    mIPSI = histcounts(outputIPSI, outputIT);
    FR_ipsi(iNeuron,:) = mIPSI/cfg.dt/length(ipsi_saccades);
    
    %% Calculate CONTRA Saccades
    [outputCONTRA, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, contra_saccades, 'doPlot', doPlot);
    mCONTRA = histcounts(outputCONTRA, outputIT);
    FR_contra(iNeuron,:) = mCONTRA/cfg.dt/length(contra_saccades);
    
    popdir;
end
A = outputIT;
binCenters = (A(:,1:end-1) + A(:,2:end)) / 2;
%% IPSI Saccades
[FR_ipsi_smooth, ~] = smoothdata(FR_ipsi,2);
[mT, iT] = max(FR_ipsi,[], 2);
FR_ipsi_norm = FR_ipsi./mT;       % normalized by FR
[~, sIPSI] = sort(iT);             % sIPSI = sorted, unnormalized
[IPSI_normSmooth, ~] = smoothdata(FR_ipsi_norm,2);

%% CONTRA Saccades
[FR_contra_smooth, ~] = smoothdata(FR_contra,2);
[mN, iN] = max(FR_contra,[], 2);
FR_contra_norm = FR_contra./mN;       % normalized by FR
[~, sCONTRA] = sort(iN);             % sIPSI = sorted, unnormalized
[CONTRA_normSmooth, ~] = smoothdata(FR_contra_norm,2);

Z.FR_ipsi = FR_ipsi;
Z.FR_contra = FR_contra;
Z.FR_ipsi_smooth = FR_ipsi_smooth;
Z.FR_contra_smooth = FR_contra_smooth;

Z.FR_ipsi_norm = FR_ipsi_norm;
Z.FR_contra_norm = FR_contra_norm;

Z.IPSI_normSmooth = IPSI_normSmooth;
Z.CONTRA_normSmooth = CONTRA_normSmooth;

Z.outputIT = outputIT;
Z.binCenters = binCenters;
Z.cfg = cfg;
Z.cellname = cellname;

Z.sIPSI = sIPSI;
Z.sCONTRA = sCONTRA;

%% Plot It
if cfg.doPlot
    clf
    subplot(1,2,1)
    imagesc(Z.binCenters, 1:length(Z.IPSI_normSmooth), Z.IPSI_normSmooth(Z.sIPSI,:));
    title('IPSI Saccades')
    xlabel('time (sec)')
    ylabel('neuron #')
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'LineStyle', '-', 'Color', 'r', 'LineWidth', LineWidth)
    
    subplot(1,2,2)
    imagesc(Z.binCenters, 1:length(Z.CONTRA_normSmooth), Z.CONTRA_normSmooth(Z.sCONTRA,:));
    title('CONTRA Saccades')
    xlabel('time (sec)')
    ylabel('neuron #')
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'LineStyle', '-', 'Color', 'r', 'LineWidth', LineWidth)
end







