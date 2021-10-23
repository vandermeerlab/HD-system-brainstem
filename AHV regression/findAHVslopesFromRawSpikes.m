function [X] = findAHVslopesFromRawSpikes(varargin)

% JJS. 3/2021. Calculate the slopes of AHV tuning curves using all spikes (i.e. not binned).

% Inputs
%   S: structure cell array of ts for spike trains

% Outputs
%   slopes: Rsq values for correlation btwn FR and AHV, divided into segments (negative AHV values, positive values, all AHV values)

% cfg_def = [];
% cfg = ProcessConfig(cfg_def,cfg_in);
TextSize = 12; 
LineWidth = 3;
FontSize = 16;
doPlot = 0;
process_varargin(varargin); 

cellCounter = 0;
fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    [S] = LoadSpikesJeff;
    fc = FindFiles('*.t', 'CheckSubdirs', 0);
    % get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.subsample_factor = 10;
    [AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, S);
    AHV_dt = median(diff(AHV_tsd.tvec));
    
    % for iCell = 1:length(S.t)
    %     vq = interp1(AHV_tsd.tvec, AHV_tsd.data, S.t{1});
    %     Sint{iCell} = tsd(S.t{1}, vq); % make a tsd of spike times and interpolated AHV values for each cell
    % end
    cfg_Q = [];
    cfg_Q.smooth = 'gauss';
    cfg_Q.gausswin_sd = 0.05;
    cfg_Q.dt = AHV_dt;
    F = MakeQfromS(cfg_Q, S); % convert to FR
    % convert to FR
    F.data = F.data ./ cfg_Q.dt;
    % find FR corresponding to each AHV sample
    F_idx = nearest_idx3(AHV_tsd.tvec, F.tvec);
    AHV_F = F.data(:,F_idx);
    
    x = AHV_tsd.data'; xR = x;
    xR(:,2) = ones(length(x),1);
    for iCell = 1:length(S.t)
        cellCounter = cellCounter + 1;
        %% All AHV values
        y = AHV_F(iCell, :); y = y';
        [~,~,~,~,statsAll] = regress(y,xR);
        a = polyfit(x, y, 1);
        %     R = corrcoef(x(~idx)',y(~idx)');
        coeffsAll(cellCounter,1) = a(1); coeffsAll(cellCounter,2) = a(2);
        rsqAll(cellCounter) = statsAll(1);
        pAll(cellCounter) = statsAll(3);
        
        %% Positive AHV values
        posIndex = x >= 0;
        [~,~,~,~,statsPos] = regress(y(posIndex),xR(posIndex,:));
        b = polyfit(x(posIndex'), y(posIndex'), 1);
        %     R = corrcoef(x(~idx & posIndex'),y(~idx & posIndex'));
        coeffsPos(cellCounter,1) = b(1); coeffsPos(cellCounter,2) = b(2);
        rsqPos(cellCounter) = statsPos(1);
        pPos(cellCounter) = statsPos(3);
        
        %% Negative AHV values
        negIndex = x <= 0;
        [~,~,~,~,statsNeg] = regress(y(negIndex),xR(negIndex,:));
        c = polyfit(x(negIndex'), y(negIndex'), 1);
        %     R = corrcoef(x(~idx & negIndex'),y(~idx & negIndex'));
        coeffsNeg(cellCounter,1) = c(1); coeffsNeg(cellCounter,2) = c(2);
        rsqNeg(cellCounter) = statsNeg(1);
        pNeg(cellCounter) = statsNeg(3);
        
        if doPlot == 1
            clf
            [filepath, name, ext] = fileparts(fc{iCell});
            %% Plot with Single Line
            subplot(1,2,2)
            hold on
            plot(AHV_tsd.data, AHV_F(iCell,:), '.', 'MarkerSize', .5);
            d = axis;
            xlim([min(x) max(x)])
            h3 = lsline;
            set(h3(1),'color', 'red')
            set(h3(1),'LineWidth',LineWidth)
            axis([-200 200 d(3) d(4)])
            %             p3 = polyfit(get(h3,'xdata'),get(h3,'ydata'),1);
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(0, d(4)/2, strcat('rsq =', num2str(rsqAll(cellCounter))), 'FontSize', TextSize)
            text(0, d(4)/2-20, strcat('p =', num2str(pAll(cellCounter))), 'FontSize', TextSize)
            set(gca, 'FontSize', 16)
            
            %% Plot with separte Lines for + and - AHV
            subplot(1,2,1)
            hold on
            plot(AHV_tsd.data(posIndex'), AHV_F(iCell,posIndex'), '.', 'MarkerSize', .5);
            xlim([min(x(posIndex)) max(x(posIndex))])
            h = lsline;
            set(h(1),'color', 'red')
            set(h(1),'LineWidth',LineWidth)
            %             p1 = polyfit(get(h,'xdata'),get(h,'ydata'),1);
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(50, d(4)/2, strcat('rsq =', num2str(rsqPos(cellCounter))), 'FontSize', TextSize)
            text(50, d(4)/2-20, strcat('p =', num2str(pPos(cellCounter))), 'FontSize', TextSize)
            
            z2 = plot(AHV_tsd.data(negIndex'), AHV_F(iCell,negIndex'), '.', 'MarkerSize', .5);
            xlim([min(x(negIndex)) max(x(negIndex))])
            h2  = refline(c(1), c(2));
            set(h2(1),'color', 'red')
            set(h2(1),'LineWidth',LineWidth)
%             set(h2(1),'color', z2.Color)
            %             h2 = lsline;
            %             p2 = polyfit(get(h2,'xdata'),get(h2,'ydata'),1);
%             xlim auto
            axis([-200 200 d(3) d(4)])
            title(name)
            xlabel('AHV (deg/sec)', 'FontSize', FontSize)
            ylabel('Firing Rate (Hz)', 'FontSize', FontSize)
            text(-100, d(4)/2, strcat('rsq =', num2str(rsqNeg(cellCounter))), 'FontSize', TextSize)
            text(-100, d(4)/2-20, strcat('p =', num2str(pNeg(cellCounter))), 'FontSize', TextSize)
            set(gca, 'FontSize', 16)
            pause
        end
    end
    popdir;
end

% plot(AHV_tsd.data, AHV_F(iCell,:), '.', 'MarkerSize', .5);


X.rsqAll = rsqAll; % R squared values for linear fit of all AHV bins
X.pAll = pAll; % p value for linear fit of all AHV bins
X.coeffsAll  = coeffsAll; % slope (1st column) and y-intercept (2nd column) for all AHV bins
X.rsqPos = rsqPos; % R squared values for linear fit of positive AHV bins
X.pPos = pPos; % p value for linear fit of positive AHV bins
X.coeffsPos = coeffsPos; %  slope (1st column) and y-intercept (2nd column) for positive AHV bins
X.rsqNeg = rsqNeg; % R squared values for linear fit of negative AHV bins
X.pNeg = pNeg; % p value for linear fit of negative AHV bins
X.coeffsNeg = coeffsNeg; %  slope (1st column) and y-intercept (2nd column) for negative AHV bins
