function [spiketrain_1,spiketrain_2, clu_1_toUse, clu_2_toUse, TCs, n, S] = senzai_extract_spike_train(mouseID)
%2024-10-34. JJS. Imports the spike train data and tuning curves from the selected mouse from the data folder. 
%   Import the data. Divide by 20000 to convert to seconds. Get the cluster identities for each spike. 
%   Inputs:
%           mouseID     - this is a string that identifies the mouse name. The options are '77c','79c','83b','85b','96d','98b'
%   Outputs: 
%           spiketrain1 - spiketrain from the first shank
%           spiketrain2 - spiketrain from the second shank
%           clu_1       - cluster identities from the first shank
%           clu_1       - cluster identities from the second shank
%           TCs.Maps    - tuning curves, 1 x nNeurons cell array 
%           TCs.z_ppln  - 1 x nNeurons double. Z-scores. [not sure how this was computed, but I use it later] 
%           n           - structure with fields for total cells, and shank1/2 num cells, etc. 
%           S           - structure with S.t{1:nNeurons} spike trains  

%% Load TCs
filename = strcat('YutaTest', mouseID, '_HeadDirection_OpenField.mat'); 
disp(strcat('importing data from', filename))
TCs = importdata(filename); % data strcuture that contains Maps, kappa_ppln, pval_ppln, z_ppln.
% Maps contains the tuning curves. pval_ppln is presumably a p value for some significance test for head direction selectivity. Not sure what kappa is.
n.numCells = length(TCs.Maps); sprintf('This mouse has %d neurons.', n.numCells)

%% Load Spikes
spiketrain_1 = importdata(strcat('YutaTest',mouseID,'.res.1')); % Shank 1 spike train
spiketrain_1 = spiketrain_1/20000; % divide by the sampling rate (20kHz) to get seconds
spiketrain_2 = importdata(strcat('YutaTest',mouseID,'.res.2')); 
spiketrain_2 = spiketrain_2/20000;

%% Load Clusters
filename = strcat('YutaTest', mouseID, '.clu.1');
disp(strcat('importing data from', filename))
clu_1 = importdata(filename);  % 'clu_1' here means the clusters from shank 1. 'clu_2' means the clusters from shank 2.
clu_1_toUse = clu_1(2:end); % remove the first value. Now it should equal length(spiketrain). 
assert(length(spiketrain_1) == length(clu_1_toUse))

clu_1_unique = unique(clu_1_toUse);
temp1 = ismember(clu_1_unique, 0:1); % 0 = artifact and 1 = noise. Remove these from consideration.git 
n.shank1_numCells = sum(temp1==0); % how many single units (clusters) are there?
n.clu_1_IDsToUse = clu_1_unique(clu_1_unique ~= 0 & clu_1_unique ~= 1); % Use clusters that are not zero or one.

filename = strcat('YutaTest', mouseID, '.clu.2');
disp(strcat('importing data from', filename))
clu_2 = importdata(filename);  % 'clu_1' here means the clusters from shank 1. 'clu_2' means the clusters from shank 2.
clu_2_toUse = clu_2(2:end); % remove the first value. Now it should equal length(spiketrain). 
assert(length(spiketrain_2) == length(clu_2_toUse))

clu_2_unique = unique(clu_2_toUse);
temp2 = ismember(clu_2_unique, 0:1); % 0 = artifact and 1 = noise. Remove these from consideration.git 
n.shank2_numCells = sum(temp2==0); % how many single units (clusters) are there?
n.clu_2_IDsToUse = clu_2_unique(clu_2_unique ~= 0 & clu_2_unique ~= 1); % Use clusters that are not zero or one.

assert(n.numCells == n.shank1_numCells + n.shank2_numCells)

%% Break up concatenated Spike Train into neurons  
clear S
counterToUse = 0;
S.type = 'ts'; % ts stands for 'timestamp' data, and is a data format in Mvdm lab. 
for iC = 1: n.numCells  % combine spike trains from each file into one structure
    if iC <= n.shank1_numCells   % use spiketrain_1 for the first x neurons. 
        S.t{iC} = spiketrain_1(clu_1_toUse == n.clu_1_IDsToUse(iC));
    elseif iC > n.shank1_numCells && iC <= n.numCells  % use spiketrain_2 for the remaining neurons. 
        counterToUse = counterToUse + 1;
        S.t{iC} = spiketrain_2(clu_2_toUse == n.clu_2_IDsToUse(counterToUse));
    else
        warning('problem with neuron count')
    end
    S.label{iC} = mouseID;
end
