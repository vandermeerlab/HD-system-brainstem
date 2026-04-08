function [neuronList, X] = calc_AHV_score_IPSI_CONTRA(cfg_in, tfilelist)
% JJS. 2024-10-17.  This function takes a cell array of neuron paths (tfilelist) and calculates the best line fit (r) and slope (in spk/s/deg/s), and AHV score
%                   defined as r x slope x 100 for each turn direction (CW & CCW), and taking the larger of those 2 values. AHV values at +-5 deg/s AHV are omitted.
% Inputs:           cfg_in      - structure with configuration variables
%                   tfilelist   - cell array. each cell is the path for one neuron, like {;'C:\Jeff\U01\datatouse\M039\M039-2020-08-21-1\M039-2020-08-21-1-TT01_1.t'}

% Outputs:          neuronList  - a list of neuron IDs that were analyzed
%                   AHV_r       - 2 x n neuron double of correlation r values. Column 1 is IPSI side. Column 2 is CONTRA side.
%                   AHV_slope   - 2 x n neuron double of best fit line slopes. Units are spk/s/deg/sec. Column 1 is IPSI side. Column 2 is CONTRA side.
%                   AHV_score   - 2 x n neuron double of correlation AHV scores. As defined by Taube (personal communication). Column 1 is IPSI side. Column 2 is CONTRA side.
%                   AHV_index   - 1 x n neuron double of whichever AHV_score is greater for each neuron.

%                       rsqAll:     [1×n double]     Rsq statistic for each neuron, fit with a single line
%                       pAll:       [1×n double]     p-values for each neuron, fit with a single line
%                       coeffsAll:  [n×2 double]     1st row = slope (Beta value); 2nd row = y-intercept. Single line fit
%                       rsqPos:     [1×n double]     Rsq for positive part of tuning curve (above zero deg/sec)
%                       pPos:       [1×n double]     p-values for positive part of tuning curve
%                       coeffsPos:  [n×2 double]     coeffs for positive part of tuning curve
%                       rsqNeg:     [1×n double]     Rsq for negative part of tuning curve
%                       pNeg:       [1×n double]     p-values for negative part of tuning curve
%                       coeffsNeg:  [n×1 double]     coeffs for negative part of tuning curve

doPlot = 0;
doSave = 0;
doPause = 0;
textFontSize = 15;
MarkerSize = 15;
FontSize = 20;
sessCounter = 0;

cfg_def.doSmooth = 0;
cfg_def.medianwindow = 10;
cfg_out = ProcessConfig2(cfg_def, cfg_in);

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
        sessCounter = sessCounter + 1;
    end
    EvalKeys
    neuronID = strcat(neuron_to_use, ext);
    cfg_spikes.fc = {neuronID};  % do one neuron at a time
    myCell = LoadSpikesJeff(cfg_spikes);
    cellname{iNeuron} = neuronID; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.)
    session_tfiles = FindFiles('*.t');  % tfiles for this session only. Need to find which number this cell is within this session. 
    [pathSession, tNameSession, ext] = fileparts(session_tfiles); 
    for iCell = 1:length(tNameSession)
        if iscell(tNameSession)
            tNameToUse = tNameSession{iCell};
        else
            tNameToUse = tNameSession;
        end
        index(iCell) = strcmp(neuron_to_use, tNameToUse);
    end
    neuronOrder = find(index); assert(~isempty(neuronOrder));
    
    %% get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.minOcc = 50; % the units for this are in samples. B/c sampling is 50Hz, 10 samples = 200 ms, for instance.
    cfg_AHV.subsample_factor = 10;
    cfg_AHV.nBins = 67; % 67 bins from -200 to 200 works out to roughly 6 degree/s bins.
    cfg_AHV.binEdges = {linspace(-200, 200, cfg_AHV.nBins)};
    
    [AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, myCell); % position of 'cfg' and 'S' got reversed at some point
    AHV_dt = median(diff(AHV_tsd.tvec));
    
    if cfg_out.doSmooth
        xdata = smoothdata(tc_out.tc, "movmedian", cfg_out.medianwindow); % moving median smoothing technique with a 10 sample (60 deg/s) moving window.
    else
        xdata = tc_out;
    end
    %% Build the matrix 
    x = xdata.binCenters; x = x';
    y = xdata.tc;    y = y';
    idx = isnan(y);
    x(:,2) = ones(length(x),1);
    
    hemisphere_to_use = ExpKeys.Hemisphere{neuronOrder};
    if strcmp(hemisphere_to_use, 'L')
        IPSIindex = xdata.binCenters > 0;  %   binCenters > 0 correspond to CCW [LEFTWARD] rotation. LEFT = IPSI for a left hemisphere neuron. 
        CONTRAindex = xdata.binCenters < 0; %           these correspond to CW [RIGHTWARD] rotation. RIGHT = CONTRA for a left hemisphere neuron. 
    elseif strcmp(hemisphere_to_use, 'R')
        IPSIindex = xdata.binCenters < 0;  %   negative values = CW = RIGHTWARD = IPSI for a right hemisphere neuron  
        CONTRAindex = xdata.binCenters > 0; %  positivate values = CCW = LEFTWARD = CONTRA for a right hemisphere neuron 
    else
        error('hemisphere error in ExpKeys')
    end
    
    %% Corr for IPSILATERAL
    [bIPSI,~,~,~,statsIPSI] = regress(y(IPSIindex),x(IPSIindex,:));
    X.bIPSI(iNeuron) = bIPSI(1);
    X.rsqIPSI(iNeuron) = statsIPSI(1);
    X.pIPSI(iNeuron) = statsIPSI(3);
    temp2 = corrcoef(x(IPSIindex' & ~idx,1),y(IPSIindex' & ~idx));
    X.rIPSI(iNeuron) = temp2(1,2);
    ahv_score_IPSI(iNeuron) = bIPSI(1)* X.rIPSI(iNeuron) *100; X.ahv_score_IPSI(iNeuron) = ahv_score_IPSI(iNeuron);
    
    %% Corr for CONTRALATERAL
    [bCONTRA,~,~,~,statsCONTRA] = regress(y(CONTRAindex),x(CONTRAindex,:));
    X.bCONTRA(iNeuron) = bCONTRA(1);
    X.rsqCONTRA(iNeuron) = statsCONTRA(1);
    X.pCONTRA(iNeuron) = statsCONTRA(3);
    temp3 = corrcoef(x(CONTRAindex' & ~idx,1),y(CONTRAindex' & ~idx));
    X.rCONTRA(iNeuron) = temp3(1,2);
    ahv_score_CONTRA(iNeuron) = bCONTRA(1)* X.rCONTRA(iNeuron) *100; X.ahv_score_CONTRA(iNeuron) = ahv_score_CONTRA(iNeuron);
    
    ahv_score(iNeuron) =  max([ahv_score_IPSI(iNeuron)  ahv_score_CONTRA(iNeuron)]); X.ahv_score(iNeuron) = ahv_score(iNeuron); 
    %% Plot It
    if doPlot; clf
        % All
        subplot(1,3,1)
        plot(tc_out.binCenters(1:end), xdata.tc(1:end), '.', 'MarkerSize', MarkerSize);
        title(cellname{iNeuron})
        xlabel('AHV (deg/s)')
        ylabel('FR (Hz)')
        set(gca, 'FontSize', FontSize)
        c = axis;
        axis([-200 200 0 c(4)]);
        line([0 0], [0 c(4)], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
        text(-180, .25*c(4), strcat('ahv Score=', num2str(round(ahv_score(iNeuron),2))), 'FontSize', textFontSize);
        
        % CW
        subplot(1,3,2)
        plot(tc_out.binCenters(1:33), xdata.tc(1:33), '.', 'MarkerSize', MarkerSize);
        lsline
        axis([-200 200 0 c(4)]);
        set(gca, 'FontSize', FontSize)
        xlabel('AHV (deg/s)')
        ylabel('FR (Hz)')
        line([0 0], [0 c(4)], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
        text(-180, .25*c(4), strcat('CW score=', num2str(round( ahv_score_CW(iNeuron),2))), 'FontSize', textFontSize);
        text(-180, .20*c(4), strcat('Rsq= ', num2str(round(X.rsqCW(iNeuron),2))), 'FontSize', textFontSize);
        text(-180, .15*c(4), strcat('p= ', num2str(round(X.pCW(iNeuron),12))), 'FontSize', textFontSize);
        title(num2str(iNeuron))
        
        % CCW
        subplot(1,3,3)
        plot(tc_out.binCenters(34:66), xdata.tc(34:66), '.', 'MarkerSize', MarkerSize);
        lsline
        set(gca, 'FontSize', FontSize)
        xlabel('AHV (deg/s)')
        ylabel('FR (Hz)')
        axis([-200 200 0 c(4)]);
        line([0 0], [0 c(4)], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
        text(20, .25*c(4), strcat('CCW score=', num2str(round( ahv_score_CCW(iNeuron),2))), 'FontSize', textFontSize);
        text(20, .20*c(4), strcat('Rsq= ', num2str(round(X.rsqCCW(iNeuron),2))), 'FontSize', textFontSize);
        text(20, .15*c(4), strcat('p= ', num2str(round(X.pCCW(iNeuron),12))), 'FontSize', textFontSize);
        
        set(gcf,'Position',[50 500 2000 500])
        hfig = tightfig(gcf);
        if doPause
            disp('press any key to continue')
            pause
        end
        if doSave
            disp('saving figure')
            cd('D:\Jeff\U01\analysis\AHV score\ahv score figs')
            saveas(gcf, neuron_to_use, 'png');  % ' D:\Jeff\U01\analysis\AHV score\ahv score figs'
        end
    end
    
end















