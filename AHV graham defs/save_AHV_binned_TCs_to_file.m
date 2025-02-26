function [TC_out, TC_out_cat, normTC, normTCsmoothed, neuronList, binCenters, sortmax, normTCsorted, TC_out_occ_cat, hemisphere] = save_AHV_binned_TCs_to_file(tfilelist, varargin)
%2024-11-26. JJS.  calculates each AHV TC (Taube params) and saves to session folder

doPlotCW = 0;
doPlotIPSI = 1;
doSave = 1;
if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
startNeuron = 1;
endNeuron = length(tfilelist);
process_varargin(varargin)
TC_out_cat = [];
TC_out_occ_cat = [];

tic
for iNeuron = startNeuron : endNeuron
    [path, neuron_to_use, ext] = fileparts(tfilelist{iNeuron});
    neuronList{iNeuron} = neuron_to_use;
    disp(neuron_to_use)
    if strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder
        pushdir(path);
        %         SSN = HD_GetSSN;
    end
    neuronID = strcat(neuron_to_use, ext);
    cfg_spikes.fc = {neuronID};  % do one neuron at a time
    myCell = LoadSpikesJeff(cfg_spikes);
    
    % Identify the hemisphere for this neuron 
    fc = FindFiles('*.t');
    [~, neurons_in_this_session, ~] = fileparts(fc);
    withinSessionNeuronNumber = find(strcmp(neuron_to_use, neurons_in_this_session));
    assert(~isempty((withinSessionNeuronNumber))); 
    
    EvalKeys;
    hemisphere{iNeuron} = ExpKeys.Hemisphere{withinSessionNeuronNumber};
    %     cellname{iNeuron} = neuronID; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.)
    
    %% get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.minOcc = 100; % the units for this are in samples. Sampling rate for rotation encoder = 200Hz. 1 sample = 5ms. This is the cutoff in Jalina's paper. They accepted anything with at least 0.5 sec data.
    cfg_AHV.subsample_factor = 10;
    cfg_AHV.nBins = 69; % 6 degree/s bins from -204 to 204 works out to 6 degree/s bins. Same as Taube lab papers. From -90 to +90 is 31 bins.
    cfg_AHV.binEdges = {linspace(-204, 204, cfg_AHV.nBins)};  % cfg_AHV.binEdges(36) == 0.
    cfg_AHV.binSize = median(diff(cfg_AHV.binEdges{1}));
    
    [AHV_tsd, tc_out, ~]  = AHV_tuning(cfg_AHV, myCell); % position of 'cfg' and 'S' got reversed at some point
    cfg_AHV.occ_dt = median(diff(AHV_tsd.tvec));
    
    binCenters = tc_out.binCenters;
    TC_out{iNeuron} = tc_out;
    
    TC_out_cat = vertcat(TC_out_cat, TC_out{iNeuron}.tc);
    TC_out_occ_cat = vertcat(TC_out_occ_cat, TC_out{iNeuron}.occ_hist);
    
    popdir;
end
normTC = TC_out_cat./max(TC_out_cat,[],2);
normTCsmoothed = smoothdata(normTC);
[~, maxi] = max(normTC, [], 2);
[~, sortmax] = sort(maxi);
normTCsorted = normTC(sortmax,:);


%% Calculate an IPSI/CONTRA version according to whether the neuron was left or right hemisphere 
hemisphere = hemisphere';
for iRow = 1: endNeuron
    if strcmp(hemisphere{iRow}, 'L')
        ICmatrix(iRow,:) = TC_out_cat(iRow,:);  % IC stands for 'ipsi/contra'
    elseif strcmp(hemisphere{iRow}, 'R')
        ICmatrix(iRow,:)  = fliplr(TC_out_cat(iRow,:));
    else
        error('hemisphere designation is neither left nor right')
    end
end

normIC = ICmatrix./max(ICmatrix,[],2);
normICsmoothed = smoothdata(normIC);
[~, maxi] = max(normIC, [], 2);
[~, sortmax] = sort(maxi);
normICsorted = normIC(sortmax,:);

if doSave
    disp('saving')
    cd('C:\Jeff\U01\datatouse\NPH analyses');
    save(strcat('Funnel plot data-', neuronList{iNeuron}));
    popdir;
end
toc
if doPlotCW
    clf
    h = imagesc(binCenters, 1:size(TC_out_cat,1), normTCsorted);
    %     h = imagesc(binCenters, 1:size(TC_out_cat,1), normTCsmoothed(sortmax,:));
    set(h, 'AlphaData', ~isnan(normTC(sortmax,:)))
    ylabel('neuron #')
    xlabel('AHV')
    set(gca, 'FontSize', 24)
    set(gca,'Color','k')
    text(-125, 100, 'CW', 'FontSize', 25, 'Color', 'w')
    text(50, 100, 'CCW', 'FontSize', 25, 'Color', 'w')
    
end
if doPlotIPSI
    clf
    h = imagesc(binCenters, 1:size(TC_out_cat,1), normICsorted);
    %     h = imagesc(binCenters, 1:size(TC_out_cat,1), normTCsmoothed(sortmax,:));
    set(h, 'AlphaData', ~isnan(normTC(sortmax,:)))
    ylabel('neuron #')
    xlabel('AHV')
    set(gca, 'FontSize', 24)
    set(gca,'Color','k')
    text(-125, 100, 'IPSI', 'FontSize', 25, 'Color', 'w')
    text(50, 100, 'CONTRA', 'FontSize', 25, 'Color', 'w')
end








%%


%     h = imagesc(binCenters, 1:size(TC_out_cat,1), normTCsmoothed(sortmax,:));
%     set(h, 'AlphaData', ~isnan(normTC(sortmax,:)))
% h = imagesc(binCenters, 1:size(ICmatrix,1), ICmatrix);

