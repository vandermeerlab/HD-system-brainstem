function [m,stats,cfg_out] = typeI_II_singleSession(cfg_in)
% JJS. 2024-06-24.
%typeI_II_singleSession.m  This function identifies for each neuron which saccade direction is IPSI and which is CONTRA. That is, it labels the 
% Side->saccade type mapping. Additionally, it examines a 1 sec. window pre-saccade and analyzes significance for phasic firing rate changes, 
% relative to a random, circularly shifted distribution of spike trains. 
%   Note: cells with significant responses can be classified as Type I or Type II. ***Need to look up which is which (as in, which is IPSI, or CONTRA)
%           Or, type II/II corresponds to slow-phase responding vs. fast-phase responing?
%           Discharge of Frontal Eye Field Neurons During Saccadic and Following Eye Movements in Unanesthetized Monkeys. Emilio Bizzi. Exp. Brain Res. 1968.
% ***I think? this is the rule for determining ipisilateral and contralateral saccades for a given neuron. NEED TO GET CONFIRMATION ON THIS 
%                   Lhem -- Leftward saccade = temporal saccade = same side? = IPSI             [L - T - I]
%                   Lhem -- Rightward saccade = nasal saccade = different side = CONTRA         [L - N - C]
%                   Rhem -- Rightward saccade =  nasal saccade = same side = IPSI               [R - N - I]
%                   Rhem -- Leftward saccade = temporal saccade = different side = CONTRA       [R - T - C]
% 
% Inputs:
%           sd - session data init. structure. Load this by typing sd = LoadSessionData([]); 
%                   sd.S = structure with fields, including .t (nSpikes x 1 double of timestamps)
%           cfg_in - config file. This can include ______
%                   cfg_in.window - response window (usually 1 sec. pre-saccade)
% Relevant keys field:
%           ExpKeys.Hemishpere = {'R'; 'L'}; % this array should be the same length as the number of neurons, and should specify 'L' or 'R' for each neuron. 
%                   If totally unsure, leave the whole array empty or include a 'NaN' value for a particular neuron that is unclear.
%           ExpKeys.HemisphereConfidence = 1;   % based on observations during recording and histology, how confident are we about which hemisphere each 
%                   TT was in? 1 = very confident. 2 = somewhat confident. 3 = not confident. Empty = no information. Note ***This should be nCells x 1 value in brackets?
% Outputs: 
%           m   
%                  .nasal     = list of whether nasal saccades are IPSI or CONTRA. Length is 1 x n neurons vector of digits. 0 = IPSI. 1 = CONTRA.
%                  .temporal  = list of whether temporal saccades are IPSI or CONTRA. Length is 1 x n neurons vector of digits. 0 = IPSI. 1 = CONTRA.
%           stats 
%                   .per = percentile values for each neuron 
%                   .sig = 0 indicates no sig. change in FR. 1 indicates sig. incr. in FR. -1 indicates sig. decrease in firing rate.
%           type -  1 x n neurons vector. For each neuron, is it an "IPSI-preferring" neuron, a "CONTRA-preferring" neuron, or neither. 
%                   -1 = IPSI sig. +1 = CONTRA sig. 0 = not sig.
 
EvalKeys;   % load the session info in ExpKeys structure
cfg_def = [];
cfg_out = ProcessConfig2(cfg_def, cfg_in);

cellnum = length(ExpKeys.Hemisphere);
stats = [];
for iCell = 1:cellnum
    if strcmp(ExpKeys.Hemisphere{iCell}, 'L')   
        m.nasal{iCell} = 'c'; % For a L hem. neuron, a nasal saccade (to the right) is CONTRA
        m.temporal{iCell} = 'i'; 
    elseif strcmp(ExpKeys.Hemisphere{iCell}, 'R')
        m.nasal{iCell} = 'i'; % For a R hem. neuron, a nasal saccade (to the right) is IPSI
        m.temporal{iCell} = 'c'; 
    else 
        error('hemisphere input is not recognized')  % ExpKeys.Hemiphere must be a cell array with values of 'R' or 'L' 
    end
end

