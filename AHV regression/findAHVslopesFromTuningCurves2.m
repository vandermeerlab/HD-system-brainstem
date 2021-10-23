function [X] = findAHVslopesFromTuningCurves2(varargin)
% JJS. 3/2021. Calculate the slopes of AHV tuning curves.
doPlot = 0;
FontSize = 14;
process_varargin(varargin);


cellCounter = 0;

fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    
    [S] = LoadSpikesJeff;
    % get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.subsample_factor = 10;
    [AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, S); % position of 'cfg' and 'S' got reversed at some point 
    AHV_dt = median(diff(AHV_tsd.tvec));
    
    fc = FindFiles('*.t', 'CheckSubdirs', 0);
    
    for iCell = 1:length(fc)
        %% Corr for Entire Tuning Curve
        cellCounter = cellCounter +1;
        x = tc_out.usr.binCenters; x = x';
        y = tc_out.tc(iCell,:);    y = y';
        x(:,2) = ones(length(x),1);
        [b,~,~,~,statsAll] = regress(y,x);
        X.b(cellCounter) = b(1);
        X.rsqAll(cellCounter) = statsAll(1);
        X.pAll(cellCounter) = statsAll(3);
        
        % Corr for Positive AHV values
        posIndex = tc_out.usr.binCenters > 0;
        [bPos,~,~,~,statsPos] = regress(y(posIndex),x(posIndex,:));
        X.bPos(cellCounter) = bPos(1);
        X.rsqPos(cellCounter) = statsPos(1);
        X.pPos(cellCounter) = statsPos(3);
        
        %% Corr for Negative AHV values
        negIndex = tc_out.usr.binCenters < 0;
        [bNeg,~,~,~,statsNeg] = regress(y(negIndex),x(negIndex,:));
        X.bNeg(cellCounter) = bNeg(1);
        X.rsqNeg(cellCounter) = statsNeg(1);
        X.pNeg(cellCounter) = statsNeg(3);
        
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
            text(50, 65, strcat('rsq =', num2str(X.rsqPos(cellCounter))))
            text(50, 45, strcat('p =', num2str(X.pPos(cellCounter))))            
            
            plot(x(negIndex,1), y(negIndex), '.');
            h2 = lsline;
            p2 = polyfit(get(h,'xdata'),get(h,'ydata'),1);
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(-100, 65, strcat('rsq =', num2str(X.rsqPos(cellCounter))))
            text(-100, 45, strcat('p =', num2str(X.pPos(cellCounter)))) 
                        
            plot(tc_out.usr.binCenters, tc_out.tc(iCell,:), '.');
            lsline
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(0, 65, strcat('rsq =', num2str(X.rsqAll(cellCounter))))
            text(0, 45, strcat('p =', num2str(X.pAll(cellCounter))))
            pause
        end
    end
    popdir;
end