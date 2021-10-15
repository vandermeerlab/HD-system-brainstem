function [z] = saccade_zscore(meanX, cellFR, varargin)

%Calculate the zscore for each time bin of a single neuron's PETH, as compared to the firing rate matrix of shuffled spike times (meanX)
%Inputs: 
%   meanX :  
%   cellFR:   1 x nBins vector of Firing Rates (PETH)
  

process_varargin(varargin);


z = (cellFR - mean(meanX))./std(meanX);
