function plot_saccade_peth_ipsi_contra(cfg_in, sd)
%2025-03-09. JJS.  Plot basic saccade peth with lines, in ipsi/contra space.

FontSize = 14;
LineWidth = 5;
colorone = [0.4660 0.6740 0.1880]; 
colortwo = 'r';

cfg_def.doPlot = 0;
cfg_def.FontSize = 16;
cfg_def.window = [-.2 .2];
cfg_def.dt = 0.01;
cfg_def.excessBounds = 1;

cfg = ProcessConfig2(cfg_def, cfg_in);


[spkTimes_N, spkTimes_T, FR_N, FR_T, bins]  = saccadePETH(cfg_def, sd);

EvalKeys;
for iNeuron = 1:length(sd.S.t)
    if strcmp(ExpKeys.Hemisphere{iNeuron}, 'R')
        ipsi{iNeuron} = FR_N{iNeuron};
        contra{iNeuron} = FR_T{iNeuron};
    elseif strcmp(ExpKeys.Hemisphere{iNeuron}, 'L')
        contra{iNeuron} = FR_N{iNeuron};
        ipsi{iNeuron} = FR_T{iNeuron};
    else
        error('issue with neuron hemisphere')
    end
    
    clf
    plot(bins(1:end-1), ipsi{iNeuron}, 'Color', colorone, 'LineWidth', LineWidth); hold on
    plot(bins(1:end-1), contra{iNeuron}, 'Color', colortwo, 'LineWidth', LineWidth)
    set(gca, 'FontSize', FontSize)
    c = axis;
    line([0 0], [0 c(4)], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    title(sd.fn{iNeuron})
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')
    legend('IPSI', 'CONTRA')
    
    disp('press any key to continue')
    pause
    
    
    
end

