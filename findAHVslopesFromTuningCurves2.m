function [rsq, p] = findAHVslopesFromTuningCurves2(varargin)
% JJS. 3/2021. Calculate the slopes of AHV tuning curves.
doPlot = 1;
FontSize = 14;
process_varargin(varargin);

rsq = []; rsqPos = []; rsqNeg = []; 
p = []; pPos = []; pNeg = []; 
cellCounter = 0;

fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    
    [S] = LoadSpikesJeff;
    % get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.subsample_factor = 10;
    [AHV_tsd, tc_out] = AHV_tuning(S, cfg_AHV);
    AHV_dt = median(diff(AHV_tsd.tvec));
    
    fc = FindFiles('*.t', 'CheckSubdirs', 0);
    
    for iCell = 1:length(fc)
        %% Corr for Entire Tuning Curve
        cellCounter = cellCounter +1;
        x = tc_out.usr.binCenters; x = x';
        y = tc_out.tc(iCell,:);    y = y';
        x(:,2) = ones(length(x),1);
        [b,~,~,~,statsAll] = regress(y,x);
        rsqAll(cellCounter) = statsAll(1);
        pAll(cellCounter) = statsAll(3);
        
        % Corr for Positive AHV values
        posIndex = tc_out.usr.binCenters > 0;
        [~,~,~,~,statsPos] = regress(y(posIndex),x(posIndex,:));
        rsqPos(cellCounter) = statsPos(1);
        pPos(cellCounter) = statsPos(3);
        
        %% Corr for Negative AHV values
        negIndex = tc_out.usr.binCenters < 0;
        [~,~,~,~,statsNeg] = regress(y(negIndex),x(negIndex,:));
        rsqNeg(cellCounter) = statsNeg(1);
        pNeg(cellCounter) = statsNeg(3);
        
        if doPlot == 1
            [filepath, name, ext] = fileparts(fc{iCell});
            subplot(1,3,1)
            hold on
            plot(x(posIndex,1), y(posIndex), '.');
            h = lsline;
            p1 = polyfit(get(h,'xdata'),get(h,'ydata'),1);
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(50, 65, strcat('rsq =', num2str(rsqPos(cellCounter))))
            text(50, 45, strcat('p =', num2str(pPos(cellCounter))))            
            
            plot(x(negIndex,1), y(negIndex), '.');
            h2 = lsline;
            p2 = polyfit(get(h,'xdata'),get(h,'ydata'),1);
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(-100, 65, strcat('rsq =', num2str(rsqPos(cellCounter))))
            text(-100, 45, strcat('p =', num2str(pPos(cellCounter)))) 
                        
            plot(tc_out.usr.binCenters, tc_out.tc(iCell,:), '.');
            lsline
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(0, 65, strcat('rsq =', num2str(rsqAll(cellCounter))))
            text(0, 45, strcat('p =', num2str(pAll(cellCounter))))
            pause
        end
    end
    popdir;
end