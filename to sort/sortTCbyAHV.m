function [isort] = sortTCbyAHV(TCs, varargin)

% 4/2021. JJS. 
% Sort tuning curves by the order of their AHV regression values and plot as imagesc. 

% INPUTS:
%
% TCs: structure containing tuning curves for all neurons. Should have fields .tc and .bins. Generated from get_TuningCurves.m. 
%       Input can be either [TC_all] (not normalized) or [TCs] (normalized). 

% OUTPUTS:
%
% sortedTCs: a structure containing sorted tuning curves for all neurons ... 
process_varargin(varargin); 

[M, I] = max(TCs.tc(:,:),[],2);

[B, isort] = sort(I);

% pcolor(TCs.bins(1,:), 1:length(TCs.tc), TCs.tc(isort, :));
imagesc(TCs.bins(1,:), 1:length(TCs.tc), TCs.tc(isort, :));

cmap = parula;
cmap(1,:) = 0;
colormap(cmap);