function [TC_out, neuronList] = save_AHV_binned_TCs_to_file(tfilelist, varargin)
%2024-11-26. JJS.  calculates each AHV TC (Taube params) and saves to session folder

if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
startNeuron = 1; 
endNeuron = length(tfilelist);
process_varargin(varargin)

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
%     cellname{iNeuron} = neuronID; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.)

    %% get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.minOcc = 100; % the units for this are in samples. Sampling rate for rotation encoder = 200Hz. 1 sample = 5ms. This is the cutoff in Jalina's paper. They accepted anything with at least 0.5 sec data.
    cfg_AHV.subsample_factor = 10;
    cfg_AHV.nBins = 69; % 6 degree/s bins from -204 to 204 works out to 6 degree/s bins. Same as Taube lab papers. From -90 to +90 is 31 bins.
    cfg_AHV.binEdges = {linspace(-204, 204, cfg_AHV.nBins)};  % cfg_AHV.binEdges(36) == 0.
    cfg_AHV.binSize = median(diff(cfg_AHV.binEdges{1}));
    
    [~, tc_out, cfg_out] = AHV_tuning(cfg_AHV, myCell); % position of 'cfg' and 'S' got reversed at some point
    
    disp('saving')
    save(strcat('AHV_binned_tuning_curve-', neuronList{iNeuron}), 'tc_out', 'cfg_out');
    
    TC_out{iNeuron} = tc_out;
    popdir;
end
toc

