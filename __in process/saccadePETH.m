function [spkTimes_N, spkTimes_T, FR_N, FR_T, bins]  = saccadePETH(cfg_in, sd)
% JJS. 2023-04-27.
% This function pulls out the spike times in a window around temporal and nasal saccades
%          Inputs:
%           - cfg_in    [struct]: contains configuration paramters
%           - S         [TS]      Spike timestamp data
%           - t         [n x T]   timestamps for events
%          Outputs:
%           - outputS   [TS] spike times, relative to t
%           - FR        1 x nBins  array of Firing Rate values
%           - bins      bins for the Firing Rate values

doPlot = 1;
cfg_def.FontSize = 16;
cfg_def.window = [-.2 .4];
cfg_def.dt = 0.01;
cfg_def.excessBounds = 1;
cfg = ProcessConfig2(cfg_def, cfg_in);

T = sd.temporalSaccades;
N = sd.nasalSaccades;

nT = length(T);
nN = length(N);
bins = linspace(cfg.window(1), cfg.window(2), diff(cfg.window)/cfg.dt+1);

for iCell = 1:length(sd.S.t)
    %% Temporal
    outputS = [];
    outputT = [];
    Stouse.t{1} = sd.S.t{iCell};
    Stouse.cfg = sd.S.cfg;
    for iT = 1:nT
        S0 = restrict(Stouse, T(iT) + cfg.window(1) - cfg.excessBounds, T(iT) + cfg.window(2) + cfg.excessBounds);
        S0 = restrict(S0, T(iT) + cfg.window(1), T(iT) + cfg.window(2));
        if length(S0.t{1}) > 0 %#ok<*ISMT>
            outputT = [outputT; repmat(iT, length(S0.t{1}),1)];
            temp = outputS;
            outputS = cat(1, temp, S0.t{1} - T(iT));   % JJS hack. Error with cat dim not consistent.
        end
    end
    [m, ~] = histcounts(outputS, bins);
    FR_T{iCell} = m / cfg.dt / length(T);
    spkTimes_T{iCell} = outputS;
    eventindexT{iCell} = outputT; 
    
    %% Nasal
    outputS = [];
    outputT = [];
    Stouse.t{1} = sd.S.t{iCell};
    for iT = 1:nN
        S0 = restrict(Stouse, N(iT) + cfg.window(1) - cfg.excessBounds, N(iT) + cfg.window(2) + cfg.excessBounds);
        S0 = restrict(S0, N(iT) + cfg.window(1), N(iT) + cfg.window(2));
        if length(S0.t{1}) > 0 %#ok<*ISMT>
            outputT = [outputT; repmat(iT, length(S0.t{1}),1)];
            temp = outputS;
            outputS = cat(1, temp, S0.t{1} - N(iT));   % JJS hack. Error with cat dim not consistent.
        end
    end
    [m, ~] = histcounts(outputS, bins);
    FR_N{iCell} = m / cfg.dt / length(N);
    spkTimes_N{iCell} = outputS;
    eventindexN{iCell} = outputT; 
    
    %% check if there are any spikes
    if isempty(outputT)
        error('No spikes')
    end
    % display
    if doPlot == 1
        clf
        %% Temporal
        subplot(2,2,1);
        title('Temporal Saccades')
        plot(spkTimes_T{iCell}, eventindexT{iCell} + 0.5, 'k.', 'MarkerSize', 5);
        xlabel('peri-event (sec)');
        ylabel('Event #');
        ylim([1 nT])
        xlim(cfg.window);
        set(gca, 'FontSize', cfg.FontSize)
        
        % bar graph
        subplot(2,2,3);
        m = histc(spkTimes_T{iCell}, bins);
        bar(bins, m / cfg.dt / length(T));
        set(gca, 'XLim', cfg.window);
        ylabel('FR (Hz)')
        xlabel('peri-event (sec)');
        set(gca, 'FontSize', cfg.FontSize)
        
        
        %% Nasal
        subplot(2,2,2);
        title('Nasal Saccades')
        plot(spkTimes_N{iCell}, eventindexN{iCell} + 0.5, 'k.', 'MarkerSize', 5);
        xlabel('peri-event (sec)');
        ylabel('Event #');
        ylim([1 nN])
        xlim(cfg.window);
        set(gca, 'FontSize', cfg.FontSize)
        
        % bar graph
        subplot(2,2,4);
        m = histc(spkTimes_N{iCell}, bins);
        bar(bins, m / cfg.dt / length(N));
        set(gca, 'XLim', cfg.window);
        ylabel('FR (Hz)')
        xlabel('peri-event (sec)');
        set(gca, 'FontSize', cfg.FontSize)
        
        disp('press any key')
        pause
    end
end
