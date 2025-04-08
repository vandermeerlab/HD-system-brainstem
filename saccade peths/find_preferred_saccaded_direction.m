function [FR_temporal_smooth, FR_nasal_smooth, temporal_normSmooth, nasal_normSmooth, sTemporal, sNasal, binCenters, cellname, cfg] = find_preferred_saccaded_direction(tfilelist, cfg_in)
% JJS. 2025-04-25. This function determines which saccade type has a higher firing rate
% Inputs:
%           tfilelist - cell array list of neurons to use
%
%
% FontSize = 30;
% LineWidth = 3;
% baseline_end = -
cfg_in.bin1 = -.05; 
cfg_in.bin2 = 0;
doPlot = 0;
doPlotMeans = 0;
cfg_def.doPlot = 0;
cfg_def.FRwindow = [-.2 .2];
cfg_def.dt = 0.005;
cfg = ProcessConfig(cfg_def,cfg_in);
cfg.pethBins = linspace(cfg.FRwindow(1), cfg.FRwindow(2), diff(cfg.FRwindow)/cfg.dt+1);
sessCounter = 0;
cd('C:\Jeff\U01\datatouse');
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
    
    %% Calculate IPSI Saccades
    cfg_in = [];
    cfg_in.doPlot = 0;
    cfg_in.doBar = 0;
    cfg_in.window = cfg.FRwindow;
    cfg_in.dt = cfg.dt;
    
    % temporal saccades
    [outputTemporal, ~, ~, outputIT, ~] = SpikePETHvdm(cfg_in, myCell, temporalSaccades, 'doPlot', doPlot);
    assert(isequal(cfg.pethBins, outputIT))
    mTemporal = histcounts(outputTemporal, outputIT);
    FR_temporal(iNeuron,:) = mTemporal/cfg.dt/length(temporalSaccades);
    
    % nasal saccades
    [outputNasal, ~, ~, outputIT, ~] = SpikePETHvdm(cfg_in, myCell, nasalSaccades, 'doPlot', doPlot);
    mNasal = histcounts(outputNasal, outputIT);
    FR_nasal(iNeuron,:) = mNasal/cfg.dt/length(nasalSaccades);
    
    popdir;
end
A = outputIT;
binCenters = (A(:,1:end-1) + A(:,2:end)) / 2;
%% temporal Saccades
[FR_temporal_smooth, ~] = smoothdata(FR_temporal,2);
[mT, iT] = max(FR_temporal,[], 2);
FR_temporal_norm = FR_temporal./mT;       % normalized by FR
[~, sTemporal] = sort(iT);             % sIPSI = sorted, unnormalized
[temporal_normSmooth, ~] = smoothdata(FR_temporal_norm,2);

%% nasal Saccades
[FR_nasal_smooth, ~] = smoothdata(FR_nasal,2);
[mN, iN] = max(FR_nasal,[], 2);
FR_nasal_norm = FR_nasal./mN;       % normalized by FR
[~, sNasal] = sort(iN);             % sIPSI = sorted, unnormalized
[nasal_normSmooth, ~] = smoothdata(FR_nasal_norm,2);

% subtraction_matrix = FR_nasal_smooth - FR_temporal_smooth;
cfg.burst_bins = cfg.bin1 <= binCenters & cfg.bin2 <= .02;
cfg.baseline_bins = binCenters < cfg.bin1;

%% This version has baseline subtraction, and takes the ***********MAX********** instead of the mean
for iNeuron = 1:size(FR_nasal_smooth,1)
    nasal_window_averaege(iNeuron) = mean(FR_nasal_smooth(iNeuron, cfg.burst_bins));
    nasal_baseline(iNeuron) = mean(FR_nasal_smooth(iNeuron, cfg.baseline_bins));
    
    temporal_window_averaege(iNeuron) = max(FR_temporal_smooth(iNeuron, cfg.burst_bins));
    temporal_baseline(iNeuron) = mean(FR_temporal_smooth(iNeuron, cfg.baseline_bins));

    % Subtract the baseline period
    nasal_subtraction(iNeuron) = nasal_window_averaege(iNeuron) - nasal_baseline(iNeuron);
    temporal_subtraction(iNeuron) = temporal_window_averaege(iNeuron) - temporal_baseline(iNeuron);
    
    subtraction_subtraction(iNeuron) = nasal_subtraction(iNeuron) - temporal_subtraction(iNeuron);
    if subtraction_subtraction(iNeuron) > 0  % nasal-preferring neuron 
        preferred_type_subtraction(iNeuron) = 1;
        type_subtraction{iNeuron} = 'nasal';
    elseif subtraction_subtraction(iNeuron) < 0 % temporal-preffering neuron
        preferred_type_subtraction(iNeuron) = 0; 
         type_subtraction{iNeuron} = 'temporal';
    else
        error('something went wrong')
    end
end


numNasal = sum(preferred_type_subtraction==1);
numTemporal = sum(preferred_type_subtraction==0);
disp(strcat('num nasal-preferring =', num2str(numNasal)))
disp(strcat('num temporal-preferring =', num2str(numTemporal)))
assert(numNasal + numTemporal == length(tfilelist));


%%
if doPlotMeans
    %%
    for iNeuron = 1:length(tfilelist)
        clf
        plot(binCenters, FR_nasal_smooth(iNeuron,:), 'LineWidth', 10); hold on
        plot(binCenters, FR_temporal_smooth(iNeuron,:), 'LineWidth', 10)
        legend('nasal', 'temporal')
        title(cellname{iNeuron})
        xlabel('time (sec)')
        ylabel('FR')
        set(gca, 'FontSize', 24)
        y = ylim;
        c = axis;
        axis([c(1) c(2) 0 c(4)])
        text(0.05, y(2)/3, type_subtraction{iNeuron}, 'FontSize', 40)
        disp(cellname{iNeuron})
        disp(strcat('nasal dFR =', num2str(round(nasal_subtraction(iNeuron),1))))
        disp(strcat('temporal dFR =', num2str(round(temporal_subtraction(iNeuron),1))))
        pause
    end
    %%
end


% %% This version just takes the mean of nasal (window) - mean of temporal (window). There is no baseline subtraction. 
% for iNeuron = 1:size(FR_nasal_smooth,1)
%     nasal_ave(iNeuron) = mean(FR_nasal_smooth(iNeuron, cfg.burst_bins));
%     temporal_ave(iNeuron) = mean(FR_temporal_smooth(iNeuron, cfg.burst_bins));
%     subtraction_score(iNeuron) = nasal_ave(iNeuron) - temporal_ave(iNeuron);
%     if subtraction_score(iNeuron) > 0  % nasal-preferring neuron 
%         preferred_type(iNeuron) = 1;
%         type{iNeuron} = 'nasal';
%     elseif subtraction_score(iNeuron) < 0 % temporal-preffering neuron
%         preferred_type(iNeuron) = 0; 
%          type{iNeuron} = 'temporal';
%     else
%         error('something went wrong')
%     end
% end










% if strcmp(ExpKeys.Hemisphere{neuronIndex}, 'R')
%     ipsi_saccade_side = {'n'}; % nasal
%     contra_saccade_side = {'t'};
% elseif strcmp(ExpKeys.Hemisphere{neuronIndex}, 'L')
%     contra_saccade_side = {'n'};
%     ipsi_saccade_side = {'t'};
% else
%     error('issue with neuron hemisphere')
% end
