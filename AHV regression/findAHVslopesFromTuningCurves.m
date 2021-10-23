function [X] = findAHVslopesFromTuningCurves(varargin)
% JJS. 3/2021. Calculate the slopes of AHV tuning curves. 
doPlot = 1;
FontSize = 16;
MarkerSize = 14;
process_varargin(varargin);
cellCounter = 0;

fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    
    [S] = LoadSpikesJeff;
    % get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.subsample_factor = 10;
    [AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, S);
    AHV_dt = median(diff(AHV_tsd.tvec));
    
    fc = FindFiles('*.t', 'CheckSubdirs', 0);
    
    for iCell = 1:length(fc)
        %% Corr for Entire Tuning Curve
        cellCounter = cellCounter +1;
        x = tc_out.usr.binCenters; x = x'; xR = x; % use this x variable for regress.m
        y = tc_out.tc(iCell,:);    y = y';
        xR(:,2) = ones(length(x),1);
        [~,~,~,~,statsAll] = regress(y,xR);
        idx = isnan(y); idx = idx';
        a = polyfit(x(~idx), y(~idx), 1);
%         R = corrcoef(x(~idx)',y(~idx)');
        coeffsAll(cellCounter,1) = a(1); coeffsAll(cellCounter,2) = a(2); 
        rsqAll(cellCounter) = statsAll(1);
        pAll(cellCounter) = statsAll(3);
        
        % Corr for Positive AHV values  i.e.  CCW turns 
        posIndex = tc_out.usr.binCenters > 0;
        [~,~,~,~,statsPos] = regress(y(posIndex),xR(posIndex,:));
        b = polyfit(x(posIndex & ~idx), y(posIndex & ~idx), 1);
        coeffsPos(cellCounter,1) = b(1); coeffsPos(cellCounter,2) = b(2); 
        rsqPos(cellCounter) = statsPos(1);
        pPos(cellCounter) = statsPos(3);
        
        %% Corr for Negative AHV values
        negIndex = tc_out.usr.binCenters < 0;
        [~,~,~,~,statsNeg] = regress(y(negIndex),xR(negIndex,:));
        c = polyfit(x(negIndex & ~idx), y(negIndex & ~idx), 1);
        coeffsNeg(cellCounter,1) = c(1); coeffsPos(cellCounter,2) = c(2);
        rsqNeg(cellCounter) = statsNeg(1);
        pNeg(cellCounter) = statsNeg(3);
        
        if doPlot == 1
            clf
            [filepath, name, ext] = fileparts(fc{iCell});
            subplot(1,2,1)
            hold on
            plot(x(posIndex,1), y(posIndex), '.', 'MarkerSize', MarkerSize);
            xlim([min(x(posIndex & ~idx)) max(x(posIndex & ~idx))]) 
            h = lsline;
%             p1 = polyfit(get(h,'xdata'),get(h,'ydata'),1);
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(50, 65, strcat('rsq =', num2str(rsqPos(cellCounter))))
            text(50, 45, strcat('p =', num2str(pPos(cellCounter))))         
            
            z2 = plot(x(negIndex,1), y(negIndex), '.', 'MarkerSize', MarkerSize);
            xlim([min(x(negIndex & ~idx)) max(x(negIndex & ~idx))]) 
            h2  = refline(c(1), c(2));
            set(h2(1),'color', z2.Color)
%             h2 = lsline;
%             p2 = polyfit(get(h2,'xdata'),get(h2,'ydata'),1);
            xlim auto
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(-100, 65, strcat('rsq =', num2str(rsqNeg(cellCounter))))
            text(-100, 45, strcat('p =', num2str(pNeg(cellCounter))))
            set(gca, 'FontSize', 16)
                        
            subplot(1,2,2) 
            plot(tc_out.usr.binCenters, tc_out.tc(iCell,:), '.', 'MarkerSize', MarkerSize);
            d = axis; 
            xlim([min(x(~idx)) max(x(~idx))]) 
            h3 = lsline;
            axis([d(1) d(2) d(3) d(4)])
%             p3 = polyfit(get(h3,'xdata'),get(h3,'ydata'),1);
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(0, 65, strcat('rsq =', num2str(rsqAll(cellCounter))))
            text(0, 45, strcat('p =', num2str(pAll(cellCounter))))
            set(gca, 'FontSize', 16)
            pause
        end
    end
    popdir;
end
X.rsqAll = rsqAll; % R squared values for linear fit of all AHV bins
X.pAll = pAll; % p value for linear fit of all AHV bins
X.coeffsAll  = coeffsAll; % slope (1st column) and y-intercept (2nd column) for all AHV bins
X.rsqPos = rsqPos; % R squared values for linear fit of positive AHV bins
X.pPos = pPos; % p value for linear fit of positive AHV bins
X.coeffsPos = coeffsPos; %  slope (1st column) and y-intercept (2nd column) for positive AHV bins
X.rsqNeg = rsqNeg; % R squared values for linear fit of negative AHV bins
X.pNeg = pNeg; % p value for linear fit of negative AHV bins
X.coeffsNeg = coeffsNeg; %  slope (1st column) and y-intercept (2nd column) for negative AHV bins