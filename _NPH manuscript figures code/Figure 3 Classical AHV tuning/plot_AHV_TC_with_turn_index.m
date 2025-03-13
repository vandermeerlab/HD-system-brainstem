function plot_AHV_TC_with_turn_index(tfilelist, turn_index)
% 2025-03-07. JJS. Plots the AHV tuning curve for each cell along with the turn index value and a text label of whether it is symmetric, asymmetric, asymmetric-unresponsive, or nonAHV
% cfg_out = ProcessConfig2(cfg_def, cfg_in);

% turn_index 
doSmooth = 0;

for iNeuron = 1:length(tfilelist)
    [path, neuron_to_use, ext] = fileparts(tfilelist{iNeuron});
    neuronList{iNeuron} = neuron_to_use;
    disp(neuron_to_use)
    
    if strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder
        pushdir(path);
        SSN = HD_GetSSN;
        EvalKeys
        sd = LoadSessionData([]);
    end
    
    if ~exist('sd', 'var')
        sd = LoadSessionData([]);
    end
    
    t = FindFiles('*.t');
    [~, neurons_in_this_session, ~] = fileparts(t);
    match = strcmp(neuron_to_use, neurons_in_this_session);
    neuronNum = find(match);
    
    cfg_in = [];
    cfg_in.doPlot = 1;
    cfg_in.smooth = 1;
    cfg_in.nBins = 60;
    cfg_in.binEdges = {linspace(-200, 200, 101)};
    cfg_in.occ_dt = median(diff(sd.AHV.tvec));
    cfg_in.minOcc = 100;  % remember that Occ is measured in samples, not in seconds.
    
    myCell = SelectTS([], sd.S, neuronNum);
    tc_out = TuningCurves(cfg_in, myCell, sd.AHV);
    
    % Add Tuning Curve
    if doSmooth
        plot(tc_out.binCenters, smoothdata(tc_out.tc), 'LineWidth', 3, 'Color', 'k');
    else
        plot(tc_out.binCenters, tc_out.tc, 'LineWidth', 3, 'Color', 'k');
    end
    xlabel('AHV (degrees/sec)')
    ylabel('FR (Hz)')
    set(groot, 'DefaultLegendInterpreter', 'none')
    c = axis;
    axis([c(1) c(2) 0 c(4)]);
    line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k')
    set(gca, 'FontSize', 16)
    disp('press any key')
    title(sd.S.label{neuronNum});
    
    [yrange] = ylim;
    [xrange] = xlim;
    ymin = yrange(1); ymax = yrange(end);
    xmin = xrange(1); xmax = xrange(end);
    text(0.75*xmax, 0.9*ymax, num2str(round(turn_index(iNeuron),2)), 'FontSize', 30, 'Color', 'r');
    
    
    pause
    
end

