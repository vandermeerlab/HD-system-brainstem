function [M, Mfr] = saccade_peth_shuffled_spike_times(myCell, t, varargin)
% Calculate significant bins for saccade peth compared to shuffled distribution. Single cell.
% spikePETH(S, t, varargin)
%
% input:
%        myCell = ts of spike times for a single neuron
%        t = Event times of interest
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
    temp = ShuffleTS(cfg, myCell);
    a.t{1} = temp.t{1};
    a.cfg = myCell.cfg; 
%     shuffmatrix(iShuff,:) = ts_in.t{1}';
    cfg_peth = [];
    cfg_peth.dt = dt;
    cfg_peth.doPlot = 0;
    cfg_peth.window = [-.2 .2];
    [outputS, outputT, outputGau, outputIT, cfg_peth] = SpikePETHvdm(cfg_peth, myCell, t);  % calculate spike peth on each shuffled spike train
    M(iShuff,:) = histc(outputS, outputIT);  % arrange the results into a matrix with histc. These are spike counts.  
    Mfr(iShuff,:) = M(iShuff,:)/dt/length(t);   % these are firing rates 
%     bar(outputIT,m/cfg.dt/length(t));
    if doPlot == 1
        clf
%         m = histc(outputS, outputIT);
        bar(outputIT,m/cfg_peth.dt/length(t));
    end
    clear a
    clear temp
end
























% Bins = window(1):dt:window(2);
% numBins = size(Bins,2);
% % For a single shuffle, loop over saccades
% h = waitbar(0, 'please wait');
% for iRow = 1:numShuff
%     waitbar(iRow/numShuff)
%     tic
%     for iT = 1:length(t)
%         for iBin = 1:numBins-1
%             X{iRow}(iT, iBin) = sum(shuffmatrix(iRow,:) > t(iT) + Bins(iBin) & shuffmatrix(iRow,:) < t(iT) + Bins(iBin+1)).*(1/dt);
%         end        
%     end
%     meanX(iRow,:) = mean(X{iRow}); % 1 x nBins, average FR
%     toc    
% end

    
    
