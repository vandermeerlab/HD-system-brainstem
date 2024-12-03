function plotAHV_tc(X, Z)
% 2024-11-26. JJS. 
% This function plots the AHV tuning curve for each neuron, with Taube style parameters, along with text boxes that describe slope, r value, and shuffle result, 
% along with whether the neuron passed all three tests (or not). 
% 6 degree bins, -90deg/s to +90deg/s. Plotting binned (average within each bin) data, now all FR samples. 
clf
        plot(tc_out.binCenters, tc_out.tc, '.', 'MarkerSize', 25); hold on
        CCWline = X.bCCW(iNeuron)*binSize*(1:length(CCWbins)) + X.CCWyInt(iNeuron);   % y = mx + b
        plot(CCWbins, CCWline)
        
        CWline = X.bCW(iNeuron)*binSize*(1:length(CWbins)) + X.CWyInt(iNeuron);  %
        
        [~, id] = max([CWline(end) X.CWyInt(iNeuron)]);
        if id == 1  % last value of negline (i.e. at y-intersect) is greater than the true y intercept
            yspan = -(CWline(end) - X.CWyInt(iNeuron));
        elseif id == 2 % last value of negline is less than the true y intercept
            yspan = abs(CWline(end) - X.CWyInt(iNeuron));
        end
        plot(CWbins, CWline + yspan)
        
        z = ylim;
        ymin = z(1);
        ymax = z(2);
        ymid = ((ymax - ymin)/2) + ymin;
        yincrement = .05*(ymax -ymin);
        
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
        
        % CW (left hand) side of plot
        text(tclose, ymid, strcat('m= ', num2str(round(bCW(iNeuron),3))), 'FontSize', 20, 'Color', ColorToUse1);
        text(tclose, ymid - yincrement, strcat('r= ', num2str(round(rCW(iNeuron),2))), 'FontSize', 20, 'Color', ColorToUse2);
        
        % CCW (right hand) side of plot 
        text(xfar, ymid, strcat('m= ', num2str(round(bCCW(iNeuron),3))), 'FontSize', 20, 'Color', ColorToUse3);
        text(xfar, ymid - yincrement, strcat('r= ', num2str(round(rCCW(iNeuron),2))), 'FontSize', 20, 'Color', ColorToUse4);
        
        xlabel('AHV')
        ylabel('FR')
        set(gca, 'FontSize', 26)
        set(gcf, 'Position', get(0, 'Screensize'));
        title(num2str(iNeuron), 'Color', titleColor)
        
        if doPause
            disp('press any key')
            pause
        end
        
end

