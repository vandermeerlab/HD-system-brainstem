function [FR, bins]  = FR_peth(cfg_in, S, t)
% JJS. 2023-04-27.
% This function pulls out the spike times in a window around temporal and nasal saccades
%          Inputs:
%           - cfg_in [struct]: contains configuration paramters
%           - S      [TS]      Spike timestamp data
%           - t      [n x T]   timestamps for events
%          Outputs:
%           - FR     1 x nBins  array of Firing Rate values 
%           - 
%           -
doPlot = 1;
cfg_def.window = [-1 2];
cfg_def.dt = 0.01;
cfg_def.excessBounds = 1;
cfg = ProcessConfig2(cfg_def, cfg_in);

nT = length(t);
outputS = [];
outputT = [];
bins = linspace(cfg.window(1), cfg.window(2), diff(cfg.window)/cfg.dt+1);

for iT = 1:nT
    S0 = restrict(S, t(iT)+cfg.window(1)-cfg.excessBounds, t(iT)+cfg.window(2)+cfg.excessBounds);
    S0 = restrict(S0, t(iT)+cfg.window(1), t(iT)+cfg.window(2));
    if length(S0.t{1}) > 0 %#ok<*ISMT>
        outputT = [outputT; repmat(iT, length(S0.t{1}),1)];
        temp = outputS;
        outputS = cat(1, temp, S0.t{1}-t(iT));   % JJS hack. Error with cat dim not consistent.
    end
end
[m, ~] = histcounts(outputS, bins);
FR = m / cfg.dt / length(t);

%% check if there are any spikes
if isempty(outputT)
    disp('No spikes')
    return
end
%% display
if doPlot ==1
    clf
    subplot(2,1,1);
    plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
    ylabel('Event #');
    ylim([1 nT])
    xlim(cfg.window);
    hold on

    % bar graph
    subplot(2,1,2);
    m = histc(outputS, bins);
    bar(bins,m/cfg.dt/length(t));
    set(gca, 'XLim', cfg.window);
    ylabel('FR (Hz)')
    xlabel('peri-event (sec)');
end
