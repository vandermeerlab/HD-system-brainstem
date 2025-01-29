function [spkTimes_N, spkTimes_T, eventindexN, eventindexT, FR_N, FR_T, bins, cfg_out]  = saccadePETH_one_neuron(cfg_in, N, T, myCell)
% JJS. 2023-04-27.
% This function pulls out the spike times in a window around temporal and nasal saccades
%          Inputs:
%           - cfg_in    [struct]: contains configuration paramters
%           - S         [TS]      Spike timestamp data
%           - t         [n x T]   timestamps for events
%           - T         temporal saccade timestamps
%           - N         nasal saccade timestamps 
%          Outputs:
%           - outputS   [TS] spike times, relative to t
%           - FR        1 x nBins  array of Firing Rate values
%           - bins      bins for the29*3 Firing Rate values

cfg_def.FontSize = 16;
cfg_def.window = [-.2 .4];
cfg_def.dt = 0.01;
cfg_def.excessBounds = 1;
cfg_out = ProcessConfig2(cfg_def, cfg_in);

nT = length(T);
nN = length(N);
bins = linspace(cfg_out.window(1), cfg_out.window(2), diff(cfg_out.window)/cfg_out.dt+1);

    %% Temporal
    outputS = [];
    outputT = [];
    for iT = 1:nT
        S0 = restrict(myCell, T(iT) + cfg_out.window(1) - cfg_out.excessBounds, T(iT) + cfg_out.window(2) + cfg_out.excessBounds);
        S0 = restrict(S0, T(iT) + cfg_out.window(1), T(iT) + cfg_out.window(2));
        if length(S0.t{1}) > 0 %#ok<*ISMT>
            outputT = [outputT; repmat(iT, length(S0.t{1}),1)];
            temp = outputS;
            outputS = cat(1, temp, S0.t{1} - T(iT));   % JJS hack. Error with cat dim not consistent.
        end
    end
    [m, ~] = histcounts(outputS, bins);
    FR_T = m / cfg_out.dt / length(T);
    spkTimes_T = outputS;
    eventindexT = outputT;
    
    %% Nasal
    outputS = [];
    outputT = [];
    for iT = 1:nN
        S0 = restrict(myCell, N(iT) + cfg_out.window(1) - cfg_out.excessBounds, N(iT) + cfg_out.window(2) + cfg_out.excessBounds);
        S0 = restrict(S0, N(iT) + cfg_out.window(1), N(iT) + cfg_out.window(2));
        if length(S0.t{1}) > 0 %#ok<*ISMT>
            outputT = [outputT; repmat(iT, length(S0.t{1}),1)];
            temp = outputS;
            outputS = cat(1, temp, S0.t{1} - N(iT));   % JJS hack. Error with cat dim not consistent.
        end
    end
    [m, ~] = histcounts(outputS, bins);
    FR_N = m / cfg_out.dt / length(N);
    spkTimes_N = outputS;
    eventindexN = outputT;
    
    %% check if there are any spikes
    if isempty(outputT)
        error('No spikes')
    end
end

% if doPlot == 1
%     clf
%     for iCell = 1:length(sd.S.t)
%         %% Temporal
%         subplot(2,2,1);
%         plot(spkTimes_T{iCell}, eventindexT{iCell} + 0.5, 'k.', 'MarkerSize', 5);
%         xlabel('peri-event (sec)');
%         ylabel('Event #');
%         ylim([1 nT])
%         xlim(cfg_out.window);
%         set(gca, 'FontSize', cfg_out.FontSize)
%         title('Temporal Saccades')
%         
%         % bar graph
%         subplot(2,2,3);
%         m = histc(spkTimes_T{iCell}, bins);
%         bar(bins, m / cfg_out.dt / length(T));
%         set(gca, 'XLim', cfg_out.window);
%         ylabel('FR (Hz)')
%         xlabel('peri-event (sec)');
%         set(gca, 'FontSize', cfg_out.FontSize)
%         
%         
%         %% Nasal
%         subplot(2,2,2);
%         plot(spkTimes_N{iCell}, eventindexN{iCell} + 0.5, 'k.', 'MarkerSize', 5);
%         xlabel('peri-event (sec)');
%         ylabel('Event #');
%         ylim([1 nN])
%         xlim(cfg_out.window);
%         set(gca, 'FontSize', cfg_out.FontSize)
%         title('Nasal Saccades')
%         
%         % bar graph
%         subplot(2,2,4);
%         m = histc(spkTimes_N{iCell}, bins);
%         bar(bins, m / cfg_out.dt / length(N));
%         set(gca, 'XLim', cfg_out.window);
%         ylabel('FR (Hz)')
%         xlabel('peri-event (sec)');
%         set(gca, 'FontSize', cfg_out.FontSize)
%         
%         disp('press any key')
%         if iCell ~= length(sd.S.t)
%             pause
%         end
%     end
% end
% 
