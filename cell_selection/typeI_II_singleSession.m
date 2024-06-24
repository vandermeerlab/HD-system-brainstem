function [m,stats] = typeI_II_singleSession(sd,cfg_in)
%typeI_II_singleSession.m  This function identifies for each neuron which saccade direction is IPSI and which is CONTRA. That is, it labels the 
% Side->saccade type mapping. Additionally, it examines a 1 sec. window pre-saccade and analyzes significance for phasic firing rate changes, 
% relative to a random, circularly shifted distribution of spike trains. 
%   Note: cells with significant responses can be classified as Type I or Type II. ***Need to look up which is which (as in, which is IPSI, or CONTRA)
%         Or, type II/II corresponds to slow-phase responding vs. fast-phase responing?
%         Discharge of Frontal Eye Field Neurons During Saccadic and Following Eye Movements in Unanesthetized Monkeys. Emilio Bizzi. Exp. Brain Res. 1968.
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
%           m - mapping between hemisphere side for each neuron and the saccade type for contra, ipsi saccades
%           stats: 
%                .per = percentile values for each neuron 
%                .sig = 0 indicates no sig. change in FR. 1 indicates sig. incr. in FR. -1 indicates sig. decrease in firing rate.
EvalKeys;   % load the session info in ExpKeys structure


cfg_out = ProcessConfig2(cfg_def, cfg_in);

cellnum = length(sd.S.t);
m = NaN(1:numcell);
stats = [];
for iCell = 1:cellnum
    




end

