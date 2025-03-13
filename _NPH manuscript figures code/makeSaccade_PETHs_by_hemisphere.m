function [Z] = makeSaccade_PETHs_by_hemisphere(cfg_in, tfilelist)
% [] = makeSaccadeHeatPlot
% JJS. 4/2021.
R_neuron_counter = 0;
L_neuron_counter = 0;

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

%% PETH params
cfg_in = [];
cfg_in.doPlot = 0;
cfg_in.doBar = 0;
cfg_in.window = cfg.FRwindow;
cfg_in.dt = cfg.dt;

if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
for iNeuron = 1:length(tfilelist)
    [path, neuron_to_use, ext] = fileparts(tfilelist{iNeuron});
    disp(neuron_to_use)
    clear ExpKeys
    if sessCounter == 0 || strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder
        pushdir(path);
%         SSN = HD_GetSSN;
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
        hemisphere = 1;
        
    elseif strcmp(ExpKeys.Hemisphere{neuronIndex}, 'L')
        contra_saccades = nasalSaccades;
        ipsi_saccades = temporalSaccades;
        hemisphere = 2;
    else
        error('issue with neuron hemisphere')
    end
    
    switch hemisphere
        case 1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT hemisphere neuron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Calculate IPSI Saccades
            R_neuron_counter = R_neuron_counter + 1;
            
            [outputIPSI, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, ipsi_saccades, 'doPlot', doPlot);
            assert(isequal(pethBins, outputIT))
            mIPSI = histcounts(outputIPSI, outputIT);
            FR_ipsi_RIGHT(R_neuron_counter,:) = mIPSI/cfg.dt/length(ipsi_saccades);
            
            %% Calculate CONTRA Saccades
            [outputCONTRA, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, contra_saccades, 'doPlot', doPlot);
            mCONTRA = histcounts(outputCONTRA, outputIT);
            FR_contra_RIGHT(R_neuron_counter,:) = mCONTRA/cfg.dt/length(contra_saccades);
            
        case 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEFT hemisphere neuron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            L_neuron_counter = L_neuron_counter + 1; 
            
            [outputIPSI, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, ipsi_saccades, 'doPlot', doPlot);
            assert(isequal(pethBins, outputIT))
            mIPSI = histcounts(outputIPSI, outputIT);
            FR_ipsi_LEFT(L_neuron_counter,:) = mIPSI/cfg.dt/length(ipsi_saccades);
            
            %% Calculate CONTRA Saccades
            [outputCONTRA, outputT, outputGau, outputIT, cfg_out] = SpikePETHvdm(cfg_in, myCell, contra_saccades, 'doPlot', doPlot);
            mCONTRA = histcounts(outputCONTRA, outputIT);
            FR_contra_LEFT(L_neuron_counter,:) = mCONTRA/cfg.dt/length(contra_saccades);

        otherwise
            error('problem with hemisphere assignment')
            
    end
    popdir;
end
A = outputIT;
binCenters = (A(:,1:end-1) + A(:,2:end)) / 2;
%% Smooth and z-score RIGHT hemisphere neurons 
% IPSI Saccades
[FR_ipsi_smooth_RIGHT, ~] = smoothdata(FR_ipsi_RIGHT,2);
[mT, iT] = max(FR_ipsi_RIGHT,[], 2);
FR_ipsi_norm_RIGHT = FR_ipsi_RIGHT./mT;       % normalized by FR
[~, sIPSI_RIGHT] = sort(iT);             % sIPSI = sorted, unnormalized
[IPSI_normSmooth_RIGHT, ~] = smoothdata(FR_ipsi_norm_RIGHT, 2);

% CONTRA Saccades
[FR_contra_smooth_RIGHT, ~] = smoothdata(FR_contra_RIGHT,2);
[mN, iN] = max(FR_contra_RIGHT,[], 2);
FR_contra_norm_RIGHT = FR_contra_RIGHT./mN;       % normalized by FR
[~, sCONTRA_RIGHT] = sort(iN);             % sIPSI = sorted, unnormalized
[CONTRA_normSmooth_RIGHT, ~] = smoothdata(FR_contra_norm_RIGHT,2);

%% Smooth and z-score LEFT hemisphere neurons 
% IPSI Saccades
[FR_ipsi_smooth_LEFT, ~] = smoothdata(FR_ipsi_LEFT,2);
[mT, iT] = max(FR_ipsi_LEFT,[], 2);
FR_ipsi_norm_LEFT = FR_ipsi_LEFT./mT;       % normalized by FR
[~, sIPSI_LEFT] = sort(iT);             % sIPSI = sorted, unnormalized
[IPSI_normSmooth_LEFT, ~] = smoothdata(FR_ipsi_norm_LEFT, 2);

% CONTRA Saccades
[FR_contra_smooth_LEFT, ~] = smoothdata(FR_contra_LEFT,2);
[mN, iN] = max(FR_contra_LEFT,[], 2);
FR_contra_norm_LEFT = FR_contra_LEFT./mN;       % normalized by FR
[~, sCONTRA_LEFT] = sort(iN);             % sIPSI = sorted, unnormalized
[CONTRA_normSmooth_LEFT, ~] = smoothdata(FR_contra_norm_LEFT,2);

% -----------------------------------------------------------------------------------------------------------------------------------------------------------------
Z.FR_ipsi_RIGHT = FR_ipsi_RIGHT;
Z.FR_ipsi_smooth_RIGHT = FR_ipsi_smooth_RIGHT; 
Z.FR_ipsi_norm_RIGHT = FR_ipsi_norm_RIGHT; 
Z.sIPSI_RIGHT = sIPSI_RIGHT; 
Z.IPSI_normSmooth_RIGHT = IPSI_normSmooth_RIGHT;

Z.FR_contra_RIGHT = FR_contra_RIGHT;
Z.FR_contra_smooth_RIGHT = FR_contra_smooth_RIGHT; 
Z.FR_contra_norm_RIGHT = FR_contra_norm_RIGHT; 
Z.sCONTRA_RIGHT = sCONTRA_RIGHT; 
Z.CONTRA_normSmooth_RIGHT = CONTRA_normSmooth_RIGHT;

Z.FR_ipsi_LEFT = FR_ipsi_LEFT;
Z.FR_ipsi_smooth_LEFT = FR_ipsi_smooth_LEFT; 
Z.FR_ipsi_norm_LEFT = FR_ipsi_norm_LEFT; 
Z.sIPSI_LEFT = sIPSI_LEFT; 
Z.IPSI_normSmooth_LEFT = IPSI_normSmooth_LEFT;

Z.FR_contra_LEFT = FR_contra_LEFT;
Z.FR_contra_smooth_LEFT = FR_contra_smooth_LEFT; 
Z.FR_contra_norm_LEFT = FR_contra_norm_LEFT; 
Z.sCONTRA_LEFT = sCONTRA_LEFT; 
Z.CONTRA_normSmooth_LEFT = CONTRA_normSmooth_LEFT;

Z.outputIT = outputIT;
Z.binCenters = binCenters;
Z.cfg = cfg;
Z.cellname = cellname;


%% Plot It
if cfg.doPlot
    clf
    subplot(2,2,1)
    imagesc(Z.binCenters, 1:size(Z.IPSI_normSmooth_LEFT,1), Z.IPSI_normSmooth_LEFT(Z.sIPSI_LEFT,:));
    title('L Hem. IPSI Saccades')
%     xlabel('time (sec)')
    ylabel('neuron #')
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'LineStyle', '-', 'Color', 'r', 'LineWidth', LineWidth)
    
    subplot(2,2,2)
    imagesc(Z.binCenters, 1:size(Z.IPSI_normSmooth_RIGHT,1), Z.IPSI_normSmooth_RIGHT(Z.sIPSI_RIGHT,:));
    title('R. Hem. IPSI Saccades')
%     xlabel('time (sec)')
    ylabel('neuron #')
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'LineStyle', '-', 'Color', 'r', 'LineWidth', LineWidth)
    
    subplot(2,2,3)
    imagesc(Z.binCenters, 1:size(Z.CONTRA_normSmooth_LEFT,1), Z.CONTRA_normSmooth_LEFT(Z.sCONTRA_LEFT,:));
    title('L Hem. CONTRA Saccades')
    xlabel('time (sec)')
    ylabel('neuron #')
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'LineStyle', '-', 'Color', 'r', 'LineWidth', LineWidth)
    
    subplot(2,2,4)
    imagesc(Z.binCenters, 1:size(Z.CONTRA_normSmooth_RIGHT,1), Z.CONTRA_normSmooth_RIGHT(Z.sCONTRA_RIGHT,:));
    title('R Hem. CONTRA Saccades')
    xlabel('time (sec)')
    ylabel('neuron #')
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [c(3) c(4)], 'LineStyle', '-', 'Color', 'r', 'LineWidth', LineWidth)
    
    
    
    
    
    
end







