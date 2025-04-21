function [Z] = makeSaccadeHeatPlotAllneurons_preferred_nonpreferred(cfg_in, type_subtraction, tfilelist)
% [] = makeSaccadeHeatPlot. JJS. 4/2021.
% Updated to be ipsi/contrat instead of temporal/nasal 3/2025.
% 2025-04-08. JJS. Changed code to make the categories preferred/non-preferred saccade direction, instead of IPS/CONTRA. 
% Relies on find_preferred_saccaded_direction.m for the preferred/non-preferred labels. 

%Inputs:    tfilelist - cell array of neuron filenames (i.e. 'C:\Jeff\U01\datatouse\M039\M039-2020-08-21-1\M039-2020-08-21-1-TT01_1.t')
%           preferred_type_subtraction - cell array of strings indicating preferred saccade direction. Either {'nasal'} or {'temporal'}

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
    
    % Designate which saccade type is preferred. This changes on a neuron by neuron basis. 
    if strcmp(type_subtraction{iNeuron}, 'nasal')
        preferred_saccades = nasalSaccades;
        nonpreferred_saccades = temporalSaccades;
    elseif strcmp(type_subtraction{iNeuron}, 'temporal')
        nonpreferred_saccades = nasalSaccades;
        preferred_saccades = temporalSaccades;
    else
        error('issue with neuron designation')
    end
    %% Calculate IPSI Saccades
    cfg_in = [];
    cfg_in.doPlot = 0;
    cfg_in.doBar = 0;
    cfg_in.window = cfg.FRwindow;
    cfg_in.dt = cfg.dt;
    
    [outputPreferred, ~, ~, outputIT_preferred, ~] = SpikePETHvdm(cfg_in, myCell, preferred_saccades, 'doPlot', doPlot);
    assert(isequal(pethBins, outputIT_preferred))
    mPreferred = histcounts(outputPreferred, outputIT_preferred);
    FR_preferred(iNeuron,:) = mPreferred/cfg.dt/length(preferred_saccades);
    
    %% Calculate CONTRA Saccades
    [outputNonPreferred, ~, ~, outputIT_nonpreferred, ~] = SpikePETHvdm(cfg_in, myCell, nonpreferred_saccades, 'doPlot', doPlot);
    mNonPreferred = histcounts(outputNonPreferred, outputIT_nonpreferred);
    FR_nonpreferred(iNeuron,:) = mNonPreferred/cfg.dt/length(nonpreferred_saccades);
    
    popdir;
end
A = outputIT_nonpreferred;
binCenters = (A(:,1:end-1) + A(:,2:end)) / 2;
%% Preferred Saccades
[FR_preferred_smooth, ~] = smoothdata(FR_preferred,2);
[mP, iP] = max(FR_preferred,[], 2);
FR_preferred_norm = FR_preferred./mP;       % normalized by max FR
[~, sPreferred] = sort(iP);             % sPreferred = sorted, unnormalized
[Preferred_normSmooth, ~] = smoothdata(FR_preferred_norm,2);

%% non-Preferred Saccades
[FR_nonpreferred_smooth, ~] = smoothdata(FR_nonpreferred,2);
[mN, iN] = max(FR_nonpreferred,[], 2);
FR_nonpreferred_norm = FR_nonpreferred./mN;       % normalized by max FR
[~, sNonPreferred] = sort(iN);             % sPreferred = sorted, unnormalized
[nonPreferred_normSmooth, ~] = smoothdata(FR_nonpreferred_norm,2);

Z.FR_preferred = FR_preferred;
Z.FR_nonpreferred = FR_nonpreferred;

Z.FR_preferred_smooth = FR_preferred_smooth;
Z.FR_nonpreferred_smooth = FR_nonpreferred_smooth;

Z.FR_preferred_norm = FR_preferred_norm;
Z.FR_nonpreferred_norm = FR_nonpreferred_norm;

Z.Preferred_normSmooth = Preferred_normSmooth;
Z.nonPreferred_normSmooth = nonPreferred_normSmooth;

Z.outputIT_nonpreferred = outputIT_nonpreferred;
Z.binCenters = binCenters;
Z.cfg = cfg;
Z.cellname = cellname;

Z.sPreferred = sPreferred;
Z.sNonPreferred = sNonPreferred;

%% Plot It
if cfg.doPlot
    clf
    subplot(1,2,1)
    imagesc(Z.binCenters, 1:length(Z.Preferred_normSmooth), Z.Preferred_normSmooth(Z.sPreferred,:));
    title('Preferred Saccades')
    xlabel('time (sec)')
    ylabel('neuron #')
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'LineStyle', '-', 'Color', 'r', 'LineWidth', LineWidth)
    
    subplot(1,2,2)
    imagesc(Z.binCenters, 1:length(Z.nonPreferred_normSmooth), Z.nonPreferred_normSmooth(Z.sNonPreferred,:));
    title('non-Preferred Saccades')
    xlabel('time (sec)')
    ylabel('neuron #')
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'LineStyle', '-', 'Color', 'r', 'LineWidth', LineWidth)
end

%% plot with matched rows, ordered by max firing bin of preferred saccade 
% figure
% subplot(1,2,1)
% imagesc(Z.binCenters, 1:length(Z.Preferred_normSmooth), Z.Preferred_normSmooth(Z.sPreferred,:));
% subplot(1,2,2)
% imagesc(Z.binCenters, 1:length(Z.nonPreferred_normSmooth), Z.nonPreferred_normSmooth(Z.sPreferred,:));

%% plot with matched rows, not normalized. 
clf
subplot(1,2,1)
imagesc(Z.binCenters, 1:length(Z.FR_preferred_smooth), Z.FR_preferred_smooth(Z.sPreferred,:));
xlabel('time (sec)')
ylabel('neuron #')
title('Preferred Saccade Direction')
set(gca, 'FontSize', 24)

subplot(1,2,2)
imagesc(Z.binCenters, 1:length(Z.FR_nonpreferred_smooth), Z.FR_nonpreferred_smooth(Z.sPreferred,:));
xlabel('time (sec)')
title('non-Preferred Saccade Direction')
set(gca, 'FontSize', 24)

a=colorbar;
a.Label.String = 'Firing Rate (Hz)';
a.Label.Position(1) = 4;
% a.PositionConstraint = 'innerposition';

% ylabel(a,'Firing Rate (Hz)','FontSize',16,'Rotation',270, 'FontSize', 24);


