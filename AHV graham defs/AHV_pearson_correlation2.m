function [X, neuronList, numSess, cfg_out] = AHV_pearson_correlation2(cfg_in, tfilelist)
%2024-11-18. JJS. Calculate the Pearson correlation for each side of the AHV tuning curve (0-90deg/s), as outlined in Graham et al., 2023. This correlation is calculated on the binned
% (averaged) tuning curve, (not on the raw data).
%  * Might make sense in future to calcaulte and store the binned tuning curves in each folder, so I don't have to recalculate each time.
% Remember: CW = negative AHV values. CCW = positive AHV values. The AHV xaxis will be CW, CCW from left to right.

% Inputs:
%       tfilelist       - a cell array with the dir of each NPH neuron. Like {C:\Jeff\U01\datatouse\M039\M039-2020-08-21-1\M039-2020-08-21-1-TT01_1.t}
%
% Outputs:
%       slopes          - the CW and CCW slopes for the best fit lines of the tuning curve (0-90deg/s). For instance,  0.025 spikes/°/s.
%                               .CW = clockwise slope. .CCW = counterclockwise slope. .CWcrit = 0 or 1 for meeting minimum slope. .CCWcrit = 0 or 1 for meetinging min slope.
%       rvalues         - Pearson r value for FR vs. AHV for the range from 0 to 90deg/s.
%                               .CW = r val. for clockwise. .CCW = r val. for counterclockwise. .CWcrit = 0 or 1 for minimum r value. .CCWcrit = 0 or 1 for min value.
%       neuronList      - easier to read list of neurons that were used in this analysis

cfg_def.min_slope = 0.025;   % from Graham et al.
cfg_def.min_r_value = 0.5;   % from Graham et al.
cfg_out = ProcessConfig2(cfg_def, cfg_in);
doPlot = 0;
doSaveFigure = 0;
doPause = 0;
numSess = 0;
tic
if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
for iNeuron = 1:length(tfilelist)
    [path, neuron_to_use, ext] = fileparts(tfilelist{iNeuron});
    neuronList{iNeuron} = neuron_to_use;
    disp(neuron_to_use)
    if strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder
        pushdir(path);
        SSN = HD_GetSSN;
        EvalKeys
        numSess = numSess + 1;
    end
    neuronID = strcat(neuron_to_use, ext);
    cfg_spikes.fc = {neuronID};  % do one neuron at a time
    myCell = LoadSpikesJeff(cfg_spikes);
    cellname{iNeuron} = neuronID; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.)
    
    %% get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.minOcc = 100; % the units for this are in samples. Sampling rate for rotation encoder = 200Hz. 1 sample = 5ms. This is the cutoff in Jalina's paper. They accepted anything with at least 0.5 sec data.
    cfg_AHV.subsample_factor = 10;
    cfg_AHV.nBins = 69; % 6 degree/s bins from -204 to 204 works out to 6 degree/s bins. Same as Taube lab papers. From -90 to +90 is 31 bins.
    cfg_AHV.binEdges = {linspace(-204, 204, cfg_AHV.nBins)};  % cfg_AHV.binEdges(36) == 0.
    cfg_AHV.binSize = median(diff(cfg_AHV.binEdges{1}));
    
    [~, tc_out, cfg_AHV] = AHV_tuning(cfg_AHV, myCell); % position of 'cfg' and 'S' got reversed at some point
    
    X.cfg_AHV = cfg_AHV;
    %     %% Determine Hemisphere
    %     EvalKeys;
    %     fc = FindFiles('*.t');
    %     [a, b, c] = fileparts(fc);
    %     d = strcat(b,c);
    %     temp = strcmp(neuronID, d);
    %     assert(sum(temp) == 1)
    %     index = find(temp);
    %     HemisphereToUse = ExpKeys.Hemisphere{index};
    
    %% Prep
    x = tc_out.binCenters; x = x';  % tc_out.binCenters(34) is the innermost negative bin [-6 to 0]. tc_out.binCenters(35) is the innermost positive bin [0 to +6].
    y = tc_out.tc;    y = y';
    idx = isnan(y);
    x(:,2) = ones(length(x),1);
    
    %% Corr for Positive AHV values (CCW)
    CCWindex = tc_out.binCenters > 0 & tc_out.binCenters <= 90;
    X.CCWbins = tc_out.binCenters(CCWindex);
    %     disp(num2str(CCWbins))
    [bCCW,~,~,~,statsCCW] = regress(y(CCWindex),x(CCWindex,:));  % b is the slope (spk/s/deg/s). stats is Rsq, F, p-value, & error variance. Rsq is the coefficient of determination.
    X.bCCW(iNeuron) = bCCW(1);
    X.CCWyInt(iNeuron) = bCCW(2);
    X.RsqCCW(iNeuron) = statsCCW(1);
    X.pCCW(iNeuron) = statsCCW(3);
    temp2 = corrcoef(x(CCWindex' & ~idx,1),y(CCWindex' & ~idx));
    X.rCCW(iNeuron) = temp2(1,2);   % this is the r value (correlation coefficient)
    
    %% Corr for Negative AHV values (CW)
    CWindex = tc_out.binCenters < 0 & tc_out.binCenters >=-90;
    X.CWbins = tc_out.binCenters(CWindex);
    %     disp(num2str(CWbins))
    [bCW,~,~,~,statsCW] = regress(y(CWindex),x(CWindex,:));
    X.bCW(iNeuron) = bCW(1);
    X.CWyInt(iNeuron) = bCW(2);
    X.RsqCW(iNeuron) = statsCW(1);
    X.pCW(iNeuron) = statsCW(3);
    temp3 = corrcoef(x(CWindex' & ~idx,1),y(CWindex' & ~idx));
    X.rCW(iNeuron) = temp3(1,2);
    
    %     if strcmp(HemisphereToUse, 'L')
    %         ipsi_shuff{iNeuron} = CW_shuff{iNeuron};
    %         contra_shuff{iNeuron} = CCW_shuff{iNeuron};
    %     elseif strcmp(HemisphereToUse, 'R')
    %         ipsi_shuff{iNeuron} = CCW_shuff{iNeuron};
    %         contra_shuff{iNeuron} = CW_shuff{iNeuron};
    %     else
    %         warning('hemisphere info not correct')
    %         ipsi_shuff{iNeuron} = {};
    %         contra_shuff{iNeuron} = {};
    %     end
    
    %% Plot it
    if doPlot
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
        text(tclose, ymid, strcat('m= ', num2str(round(X.bCW(iNeuron),3))), 'FontSize', 20, 'Color', ColorToUse1);
        text(tclose, ymid - yincrement, strcat('r= ', num2str(round(X.rCW(iNeuron),2))), 'FontSize', 20, 'Color', ColorToUse2);
        
        % CCW (right hand) side of plot
        text(xfar, ymid, strcat('m= ', num2str(round(X.bCCW(iNeuron),3))), 'FontSize', 20, 'Color', ColorToUse3);
        text(xfar, ymid - yincrement, strcat('r= ', num2str(round(X.rCCW(iNeuron),2))), 'FontSize', 20, 'Color', ColorToUse4);
        
        if abs(X.bCW(iNeuron)) >= 0.025 || abs(X.bCCW(iNeuron)) >= 0.025
            slopeIsGood(iNeuron) = 1;
        else
            slopeIsGood(iNeuron) = 0;
        end
        
        if abs(X.rCW(iNeuron)) >= 0.5 || abs(X.rCCW(iNeuron)) >= 0.5
            rIsGood(iNeuron) = 1;
        else
            rIsGood(iNeuron) = 0;
        end
        
        if slopeIsGood(iNeuron) && rIsGood(iNeuron)
            bothIsGood(iNeuron) = 1;
            titleColor = 'r';
        else
            titleColor = 'k';
            bothIsGood(iNeuron) = 0;
        end
        
        xlabel('AHV')
        ylabel('FR')
        set(gca, 'FontSize', 26)
        set(gcf, 'Position', get(0, 'Screensize'));
        title(num2str(iNeuron), 'Color', titleColor)
        
        if doPause
            disp('press any key')
            pause
        end
        
    else  % hackey. shouldn't repeat this whole code block 
        if abs(X.bCW(iNeuron)) >= 0.025 || abs(X.bCCW(iNeuron)) >= 0.025
            slopeIsGood(iNeuron) = 1;
        else
            slopeIsGood(iNeuron) = 0;
        end
        
        if abs(X.rCW(iNeuron)) >= 0.5 || abs(X.rCCW(iNeuron)) >= 0.5
            rIsGood(iNeuron) = 1;
        else
            rIsGood(iNeuron) = 0;
        end
        
        if slopeIsGood(iNeuron) && rIsGood(iNeuron)
            bothIsGood(iNeuron) = 1;
            titleColor = 'r';
        else
            titleColor = 'k';
            bothIsGood(iNeuron) = 0;
        end
    end
    %% Save it
    if doPlot && doSaveFigure
        disp('saving figure')
        cd('D:\Jeff\U01\analysis\AHV selectivity JG definition')
        saveas(gcf,strcat(num2str(iNeuron), '.png'))
    end
end
X.slopeIsGood = slopeIsGood;
X.rIsGood = rIsGood;
X.bothIsGood = bothIsGood;
toc