function [] = saveAHV_FR_scatterplot(varargin)

% 4/2021. JJS.

% This version saves the data into a single structure instead of saving it to the session folder. 

doPlot = 0;
doSave = 1;
filename = 'FRxAHV';
process_varargin(varargin)

fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    
    SSN = HD_GetSSN;
    disp(SSN)
    
    [S] = LoadSpikesJeff;
    % get AHV Tuning Curve
    cfg_AHV = [];
    cfg_AHV.subsample_factor = 10;
    [AHV_tsd, tc_out] = AHV_tuning(S, cfg_AHV);
    AHV_dt = median(diff(AHV_tsd.tvec));
    
    %% Firing Rate Scatterplot
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
    %     ymax = max(AHV_F(iCell,:));
    
    if doPlot == 1
        clf
        for iCell = 1:length(S.t)
        plot(AHV_tsd.data, AHV_F(iCell,:), '.', 'MarkerSize', .5);
        disp('waiting for user input')
        title(S.label{iCell})
        pause
        end
    end
    if doSave == 1
        save(strcat(SSN, '-', filename, '.mat'), 'AHV_tsd', 'AHV_F')
        disp('data saved')
    end
    popdir;
end