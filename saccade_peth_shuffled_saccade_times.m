function [M, Mfr] = saccade_peth_shuffled_saccade_times(myCell, t, varargin)
% Generates a matrix of shuffled spike times by randomly permuting the saccade times. 
% Calculates significant bins for saccade peth compared to shuffled distribution. Single cell.
%
% input:
%        myCell: a ts of spike times for a single neuron 
%        t:      Event times (saccade times)  
% output:
%        z score for each time bin
% parameters:
%   window = [-.2 .2]
dt = 0.01; % seconds
doPlot = 0;
numShuff = 100;
% window = [-0.2 0.2];
% dt = 0.01;
process_varargin(varargin);

cfg = [];
cfg.mode = 2;
cfg.t0 = myCell.t{1}(1);
cfg.t1 = myCell.t{1}(end);

% shuffmatrix = nan(numShuff, size(myCell.t{1},1));
% generate a matrix of shuffled spike times
h = waitbar(0, 'please wait');
for iShuff = 1:numShuff
    waitbar(iShuff/numShuff)
    disp(iShuff) 
	ISI = diff(horzcat([t(1), t]));
	tShuff = cumsum(ISI(randperm(length(ISI))));

    cfg_peth = [];
    cfg_peth.dt = dt;
    cfg_peth.doPlot = 0;
    cfg_peth.window = [-.2 .2];
    [outputS, outputT, outputGau, outputIT, cfg_peth] = SpikePETHvdm(cfg_peth, myCell, tShuff);  % calculate spike peth on each shuffled spike train
    M(iShuff,:) = histc(outputS, outputIT);  % arrange the results into a matrix with histc. These are spike counts.  
    Mfr(iShuff,:) = M(iShuff,:)/dt/length(t);   % these are firing rates 
%     bar(outputIT,m/cfg.dt/length(t));
    if doPlot == 1
        clf
%         m = histc(outputS, outputIT);
        bar(outputIT,m/cfg_peth.dt/length(t));
    end
    clear tShuff
end
% Z = (FRxBinT(1,:) - mean(Mfr(1:end-1)))/std(Mfr(1:end-1)); % calculate the z-score. Ignore last column, which is zeros