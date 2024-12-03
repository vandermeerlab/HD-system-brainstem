function [neuronID] = plotAHV_tc_ipsi_contra(X, Z, IC, tfilelist, varargin)
% 2024-11-26. JJS.
% This function uses outputs from AHV_pearson_correlation2.m 
% This function uses outputs from calc_circ_shift_AHV_IPSI_CONTRA.m 
% This function plots the AHV tuning curve for each neuron, with Taube style parameters, along with text boxes that describe slope, r value, and shuffle result,
% along with whether the neuron passed all three tests (or not).
% 6 degree bins, -90deg/s to +90deg/s. Plotting binned (average within each bin) data, now all FR samples.
doSave = 1;
if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
startNeuron = 1;
endNeuron = length(tfilelist);
doPause = 0;
process_varargin(varargin)

tic
for iNeuron = startNeuron : endNeuron
    [path, neuron_to_use, ext] = fileparts(tfilelist{iNeuron});
    neuronList{iNeuron} = neuron_to_use;
    disp(neuron_to_use)
    if strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder
        pushdir(path);
        %         SSN = HD_GetSSN;
    end
    neuronID = strcat(neuron_to_use, ext);
    filenametouse = strcat('AHV_binned_tuning_curve-', neuronList{iNeuron});
    load(filenametouse);
    binsize = median(diff(tc_out.binCenters));
    
    clf
    plot(tc_out.binCenters, tc_out.tc, '.', 'MarkerSize', 25); hold on
    ipsi_line = X.bCCW(iNeuron)*X.cfg_AHV.binSize*(1:length(X.CCWbins)) + X.CCWyInt(iNeuron);   % y = mx + b
    plot(X.CCWbins, ipsi_line)
    
    contra_line = X.bCW(iNeuron)*binsize*(1:length(X.CWbins)) + X.CWyInt(iNeuron);  %
    
    [~, id] = max([contra_line(end) X.CWyInt(iNeuron)]);
    if id == 1  % last value of negline (i.e. at y-intersect) is greater than the true y intercept
        yspan = -(contra_line(end) - X.CWyInt(iNeuron));
    elseif id == 2 % last value of negline is less than the true y intercept
        yspan = abs(contra_line(end) - X.CWyInt(iNeuron));
    end
    plot(X.CWbins, contra_line + yspan)
    
    z = ylim;
    ymin = z(1);
    ymax = z(2);
    ymid = ((ymax - ymin)/2) + ymin;
    yincrement = .05*(ymax -ymin);
    ylow = ymin + yincrement;
    
    
    t = xlim;
    xmin = t(1);
    xmax = t(2);
    tclose = 0.8*xmin;
    xfar = 0.5*xmax;
    
    if abs(X.bCW(iNeuron)) >= 0.025
        ColorToUse1 = 'r';
    else
        ColorToUse1 = 'k';
    end
    
    if abs(X.bCCW(iNeuron)) >= 0.025
        ColorToUse3 = 'r';
    else
        ColorToUse3 = 'k';
    end
    
    if abs(X.rCW(iNeuron)) >= 0.5
        ColorToUse2 = 'r';
    else
        ColorToUse2 = 'k';
    end
    
    if abs(X.rCCW(iNeuron)) >= 0.5
        ColorToUse4 = 'r';
    else
        ColorToUse4 = 'k';
    end
    
    if Z.pass(iNeuron)
        titleColor = 'r';
    else
        titleColor = 'k';
    end
    if Z.rank_CW(iNeuron) >= 95
        ColorToUse5 = 'r';
    else
        ColorToUse5 = 'k';
    end
    if Z.rank_CCW(iNeuron) >= 95
        ColorToUse6 = 'r';
    else
        ColorToUse6 = 'k';
    end 
    
    % CW (left hand) side of plot
    text(tclose, ymid, strcat('m= ', num2str(round(X.bCW(iNeuron),3))), 'FontSize', 20, 'Color', ColorToUse1);
    text(tclose, ymid - yincrement, strcat('r= ', num2str(round(X.rCW(iNeuron),2))), 'FontSize', 20, 'Color', ColorToUse2);
    text(tclose, ymid - 2*yincrement, strcat('rank=', num2str(Z.rank_CW(iNeuron)), '%'), 'FontSize', 20, 'Color', ColorToUse5)
    text(tclose, ylow, 'CW', 'FontSize', 36)
    
    % CCW (right hand) side of plot
    text(xfar, ymid, strcat('m= ', num2str(round(X.bCCW(iNeuron),3))), 'FontSize', 20, 'Color', ColorToUse3);
    text(xfar, ymid - yincrement, strcat('r= ', num2str(round(X.rCCW(iNeuron),2))), 'FontSize', 20, 'Color', ColorToUse4);
    text(xfar, ymid - 2*yincrement, strcat('rank=', num2str(Z.rank_CCW(iNeuron)), '%'), 'FontSize', 20, 'Color', ColorToUse6)
    text(xfar, ylow, 'CCW', 'FontSize', 36)
    
    xlabel('AHV')
    ylabel('FR')
    set(gca, 'FontSize', 26)
    set(gcf, 'Position', get(0, 'Screensize'));
    title(num2str(iNeuron), 'Color', titleColor)
    
    if doPause
        disp('press any key')
        pause
    end
    if doSave
        disp('saving figure')
        cd('D:\Jeff\U01\analysis\AHV selectivity JG definition\calc_circ_shift_AHV_IPSI_CONTRA\individual TC plots')
        saveas(gcf,strcat(num2str(iNeuron), '.png'))
    end
end


