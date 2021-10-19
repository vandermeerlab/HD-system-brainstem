function [M, Mfr, zT, zN] = saccade_peth_shuffled_saccade_times2(cfg_in, myCell, cellnum, temporalSaccades, nasalSaccades, varargin)
% Generates a matrix of shuffled spike times by randomly permuting the saccade times.
% Calculates significant bins for saccade peth compared to shuffled distribution. Single cell.
%
% input:
%        myCell:    a ts of spike times for a single neuron
%        cellnum:   the order of the cell in the structure S (i.e., 1, 2, 3, ...) 
%        t:         Event times (saccade times)
% output:
%        M:     matrix of shuffled spike counts
%        Mfr:   matrix of firing rates 
%        zT:    z-score for each time bin (temporal saccades)
%        zN:    z-score for each time bin (nasal saccades)

cfg_def.FRwindow = [-.2 .2];
cfg_def.dt = 0.01;
cfg = ProcessConfig2(cfg_def, cfg_in);

doPlot = 0;
numShuff = 100;
process_varargin(varargin);

[FRxBinT, FRxBinN, ~, ~, ~, ~, ~, ~, ~, binCenters, cfg_out, ~, ~] = makeSaccadeHeatPlot(cfg, temporalSaccades, nasalSaccades, 'doPlot', 0);
h = waitbar(0, 'please wait');
for iJ = 1:2  % 2 cases
    if iJ == 1; t = temporalSaccades; end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    if iJ == 2; t = nasalSaccades; end
    for iShuff = 1:numShuff
        waitbar(iShuff/numShuff)
%         disp(iShuff)
        ISI = diff(horzcat([t(1), t]));
        tShuff = cumsum(ISI(randperm(length(ISI))));
        
        cfg_peth = [];
        cfg_peth.dt = cfg.dt;
        cfg_peth.doPlot = 0;
        cfg_peth.window = cfg.FRwindow;
        [outputS, outputT, outputGau, outputIT, cfg_peth] = SpikePETHvdm(cfg_peth, myCell, tShuff);  % calculate spike peth on each shuffled spike train
        M(iShuff,:) = histc(outputS, outputIT);  % arrange the results into a matrix with histc. These are spike counts.
        Mfr(iShuff,:) = M(iShuff,:)/cfg.dt/length(t);   % these are firing rates
        %     clear tShuff
        if iJ ==1  % temporal
            zT = (FRxBinT(cellnum,:) - mean(Mfr(1:end-1)))/std(Mfr(1:end-1)); % calculate the z-score. Ignore last column, which is zeros
        elseif iJ == 2   % nasal
            zN = (FRxBinN(cellnum,:) - mean(Mfr(1:end-1)))/std(Mfr(1:end-1));
        end
    end
end

