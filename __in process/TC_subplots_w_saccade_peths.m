function TC_subplots_w_saccade_peths(sd)

%2025-01-25. JJS. Function to plot all of the relevant tuning curves for a headfixed session with eyetracking, including
%                 HD, AHV, wheel speed, saccade peths (all saccades), and saccade peths (toward/away PFD). Each TC has a corresponding
%                 occupancy plot below it. So, for HD, how much time was spend at each orientation, etc.

% Inputs:         cfg_in   -  config parameters
%                 sd       -  session data [LoadSessionData.m]

% Outputs:        HD_tc    -  head direction tuning curve
%                 AHV_tc   -  AHV tuning curve
%                 wheel_tc -  wheel speed tuning curve
Xtight = .04;
Ytight = .02;
FontSizeToUse = 12;
AHVmin = -150;
AHVmax = 150;
LineWidth = 2;

tic
%% Calculations
updownTSD = getQEupdown([]);
state_tsd = ConvertQEUpDownToState(updownTSD);
[~, wheel_tsd] = ConvertQEStatesToAngle([], state_tsd);
[d, speed, cfg_w] = ConvertWheeltoSpeed([], wheel_tsd);   % d = distance.  speed = speed of the wheel in cm/sec
speed.tvec = speed.tvec - sd.starttime;
speed.data = -speed.data; % we want forward motion to be displayed as a positive velocity

cfg_wheel = [];
cfg_wheel.doPlot = 0;
cfg_wheel.nBins = 50;
cfg_wheel.binEdges = {linspace(-5, 30, 101)};
cfg_wheel.occ_dt = median(diff(speed.tvec));
cfg_wheel.minOcc = 1;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
toc
%% Plot Each Neuron
for iCell = 1:length(sd.fc)
    clf
    %% HD TC
    % 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
    subtightplot(4,5,1,[Xtight Ytight]); hold on
    %     set(gca, 'YTick', [])
    set(gca, 'XTick', -180:30:180)
    %     grid on
    ax = gca;
    ax.XGrid = 'on';
    ax.GridLineStyle = '-';
    
    title('HD TC: mov. + stat.', 'FontSize', FontSizeToUse)
    ylabel('FR (Hz)', 'FontSize', FontSizeToUse)
    cfg_hd = [];
    cfg_hd.doPlot = 0;
    cfg_hd.smooth = 0;
    cfg_hd.nBins = 60;
    cfg_hd.binEdges = {linspace(-180, 180, 101)};
    cfg_hd.occ_dt = median(diff(sd.AHV.tvec));  % Remember that the platform encoder samples at 2kHz.
    cfg_hd.minOcc = 200;  % remember that Occ is measured in samples, not in seconds. Usually 5ms per sample, b/c the platform encoder sampling rate is 200Hz. minOcc of 200 = 1 second duration.
    
    myCell = SelectTS([], sd.S, iCell);
    tc_HD = TuningCurves(cfg_hd, myCell, sd.orientation);
    plot(tc_HD.binCenters, tc_HD.tc, 'LineWidth', 3, 'Color', 'k'); hold on
    plot(tc_HD.binCenters, smoothdata(tc_HD.tc), 'LineWidth', 3, 'Color', 'r', 'LineStyle', '--');
    
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), 'deg','HorizontalAlignment','right','VerticalAlignment','top')
    
    %% AHV TC
    % 2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
    subtightplot(4,5,2,[Xtight Ytight]);
    set(gca, 'YTick', [])
    
    cfg_AHV = [];
    cfg_AHV.doPlot = 1;
    cfg_AHV.nBins = 5;
    cfg_AHV.binEdges = {linspace(-150, 150, 50)};
    cfg_AHV.occ_dt = median(diff(sd.AHV.tvec));
    cfg_AHV.minOcc = 100;  % remember that Occ is measured in samples, not in seconds. Usually 5ms per sample, b/c the platform encoder sampling rate is 200Hz,
    
    tc_AHV = TuningCurves(cfg_AHV, myCell, sd.AHV);
    plot(tc_AHV.binCenters, tc_AHV.tc, 'LineWidth', 3, 'Color', 'k'); hold on
    plot(tc_AHV.binCenters, smoothdata(tc_AHV.tc), 'LineWidth', 3, 'Color', 'r', 'LineStyle', '--');
    c = axis;
    axis([AHVmin AHVmax 0 c(4)]);
    
    ax = gca;
    ax.XGrid = 'on';
    ax.GridLineStyle = '-';
    title('AHV TC: mov. + stat.', 'FontSize', FontSizeToUse)
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), 'deg/sec','HorizontalAlignment','right','VerticalAlignment','top')
    
    %% Wheel TC
    % 3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
    subtightplot(4,5,3,[Xtight Ytight]);
    set(gca, 'YTick', [])
    
    tc_wheel = TuningCurves(cfg_wheel, myCell, speed);
    plot(tc_wheel.binCenters, tc_wheel.tc, 'k', 'LineWidth', 3); hold on
    plot(tc_wheel.binCenters, smoothdata(tc_wheel.tc), 'LineWidth', 3, 'Color', 'r', 'LineStyle', '--');
    
    ax = gca;
    ax.XGrid = 'on';
    ax.GridLineStyle = '-';
    title('Wheel TC: mov. + stat.', 'FontSize', FontSizeToUse)
    c = axis; 
    axis([-5 c(2) 0 c(4)]);
    set(gca, 'XTick', -5:5:30)
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), 'cm/sec','HorizontalAlignment','right','VerticalAlignment','top')
    
    
    %% Pupil Position TC 
    % 4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
    subtightplot(4,5,4,[Xtight Ytight]);
    set(gca, 'YTick', [])
    
    cfg_pupil.binEdges = {linspace(-60, 60, 61)};
    cfg_pupil.occ_dt = median(diff(sd.tsdH.tvec));
    cfg_pupil.minOcc = 10;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
    tc_pupil_position = TuningCurves([], myCell, sd.tsdH);
    
    plot(tc_pupil_position.binCenters, tc_pupil_position.tc, 'LineWidth', 3, 'Color', 'k'); hold on
    plot(tc_pupil_position.binCenters, smoothdata(tc_pupil_position.tc), 'LineWidth', 3, 'Color', 'r', 'LineStyle', '--');
    c = axis;
    axis([-30 30 0 c(4)]);

    title('Pupil Position TC: mov. + stat.', 'FontSize', FontSizeToUse)
    ylabel('FR (Hz)', 'FontSize', FontSizeToUse)
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), 'pixels','HorizontalAlignment','right','VerticalAlignment','top')
    
    %% Pupil Velocity TC 
    % 555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
    subtightplot(4,5,5,[Xtight Ytight]);
    set(gca, 'YTick', [])
    cfg_pupil_velocity = [];
    cfg_pupil_velocity.binEdges = {linspace(-40, 40, 41)};
    tc_pupil_velocity.occ_dt = median(diff(sd.diffH.tvec));
    tc_pupil_velocity.minOcc = 100;  % remember that Occ is measured in samples, not in seconds. Usually 20ms per sample, b/c the camera sampling rate is 50Hz,
    
    tc_pupil_velocity = TuningCurves(cfg_pupil_velocity, myCell, sd.diffH);
    plot(tc_pupil_velocity.binCenters, tc_pupil_velocity.tc, 'LineWidth', 3, 'Color', 'k'); hold on
    plot(tc_pupil_velocity.binCenters, smoothdata(tc_pupil_velocity.tc), 'LineWidth', 3, 'Color', 'r', 'LineStyle', '--');
    title('Pupil Velocity TC: mov. + stat.', 'FontSize', FontSizeToUse)
    
    c = axis;
    axis([-30 30 0 c(4)]);
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), 'pixels/sec','HorizontalAlignment','right','VerticalAlignment','top')
    
    %% HD Occupancy
    % 666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666
    subtightplot(4,5,6,[Xtight Ytight]);
    plot(tc_HD.binCenters, tc_HD.occ_hist)
    %     set(gca, 'YTick', [])
    set(gca, 'XTick', -180:30:180)
    grid on
    title('HD Occupancy', 'FontSize', FontSizeToUse)
    ylabel('Seconds', 'FontSize', FontSizeToUse)
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), 'deg','HorizontalAlignment','right','VerticalAlignment','top')
    
    %% AHV Occupancy
    % 777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777
    subtightplot(4,5,7,[Xtight Ytight]);
    plot(tc_AHV.binCenters, tc_AHV.occ_hist)
    c = axis;
    axis([c(1) c(2) 0 15]);
    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    set(gca, 'YTick', 0:3:15)
    title('AHV Occupancy', 'FontSize', FontSizeToUse)
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), 'deg/sec','HorizontalAlignment','right','VerticalAlignment','top')
    
    %% Wheel Occupancy
    % 88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
    subtightplot(4,5,8,[Xtight Ytight]);
    plot(tc_wheel.binCenters, tc_wheel.occ_hist)
    c = axis;
    axis([-5 c(2) 0 15]);
    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    ax.GridLineStyle = '-';
    set(gca, 'YTick', 0:3:15)
    title('Wheel Occupancy', 'FontSize', FontSizeToUse)
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), 'cm/sec','HorizontalAlignment','right','VerticalAlignment','top')
    
    %% Pupil Position Occupancy 
    % 99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999
    subtightplot(4,5,9,[Xtight Ytight]);
    plot(tc_pupil_position.binCenters, tc_pupil_position.occ_hist)
    title('Pupil Position Occupancy', 'FontSize', FontSizeToUse)
    c = axis;
    axis([-30 30 0 c(4)]);
    set(gca, 'XTick', -30:5:30)
    ax = gca;
    ax.XGrid = 'on';
    ax.GridLineStyle = '-';
    ax.GridLineStyle = '-';
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), 'pixels','HorizontalAlignment','right','VerticalAlignment','top')
    
    %% Pupil Velocity Occupancy 
    % 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 
    subtightplot(4,5,10,[Xtight Ytight]);
    plot(tc_pupil_velocity.binCenters, tc_pupil_velocity.occ_hist)
    title('Pupil Velocity Occupancy', 'FontSize', FontSizeToUse)
    c = axis;
    axis([-10 10 0 30]);
    set(gca, 'XTick', -10:2:10)
    ax = gca;
    ax.XGrid = 'on';
    ax.GridLineStyle = '-';
    ax.GridLineStyle = '-';
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), 'pixels/sec','HorizontalAlignment','right','VerticalAlignment','top')

%% GET SACCADE INFO
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    [~, ~, ~, ~, nasal_timestamps_MOVING, temporal_timestamps_MOVING] = isolateManualSaccades();
    [~, ~, nasal_timestamps_STATIONARY, temporal_timestamps_STATIONARY] = isolateStationarySaccades();

    numNasal_moving = length(nasal_timestamps_MOVING);
    numTemporal_moving = length(temporal_timestamps_MOVING);
    
    numNasal_stationary = length(nasal_timestamps_STATIONARY);
    numTemporal_stationary = length(temporal_timestamps_STATIONARY);
    
    [spkTimes_N_moving, spkTimes_T_moving, eventindexN_moving, eventindexT_moving, FR_N_moving, FR_T_moving, bins_moving, cfg_out_moving]  = saccadePETH_one_neuron([], nasal_timestamps_MOVING, temporal_timestamps_MOVING, myCell);
    [spkTimes_N_stationary, spkTimes_T_stationary, eventindexN_stationary, eventindexT_stationary, FR_N_stationary, FR_T_stationary, bins_stationary, cfg_out_stationary]  = saccadePETH_one_neuron([], nasal_timestamps_STATIONARY, temporal_timestamps_STATIONARY, myCell);

    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    %% Saccade peth Nasal MOVING 
    % 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 
    subtightplot(4,5,11,[Xtight Ytight]);
    m = histc(spkTimes_N_moving, bins_moving);
    bar(bins_moving, m / cfg_out_moving.dt / numNasal_moving);
    title('Nasal Saccades: MOVING', 'FontSize', FontSizeToUse)
    ylabel('FR (Hz)', 'FontSize', FontSizeToUse); 
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), num2str(numNasal_moving),'HorizontalAlignment','right','VerticalAlignment','top')
    c = axis; line([0 0], [0 c(4)], 'LineStyle', '--', 'Color', 'r','LineWidth', 2)
    %% Saccade peth Temporal MOVING
    % 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 12 
    subtightplot(4,5,12,[Xtight Ytight]);
    n = histc(spkTimes_T_moving, bins_moving);
    bar(bins_moving, n / cfg_out_moving.dt / numTemporal_moving);
    title('Temporal Saccades: MOVING', 'FontSize', FontSizeToUse); 
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), num2str(numTemporal_moving),'HorizontalAlignment','right','VerticalAlignment','top')
    c = axis; line([0 0], [0 c(4)], 'LineStyle', '--', 'Color', 'r','LineWidth', 2)
    %% Saccade peth Nasal STATIONARY
    % 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 
    subtightplot(4,5,16,[Xtight Ytight]);
    o = histc(spkTimes_N_stationary, bins_stationary);
    bar(bins_stationary, o / cfg_out_stationary.dt / numNasal_stationary);
    title('Nasal Saccades: STATIONARY', 'FontSize', FontSizeToUse)
    ylabel('FR (Hz)', 'FontSize', FontSizeToUse); 
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), num2str(numNasal_stationary),'HorizontalAlignment','right','VerticalAlignment','top')
    c = axis; line([0 0], [0 c(4)], 'LineStyle', '--', 'Color', 'r','LineWidth', 2)
    %% Saccade peth Temporal STATIONARY
    % 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 
    subtightplot(4,5,17,[Xtight Ytight]);
    p = histc(spkTimes_T_stationary, bins_stationary);
    bar(bins_stationary, p / cfg_out_stationary.dt / numTemporal_stationary);
    title('Temporal Saccades: STATIONARY', 'FontSize', FontSizeToUse); 
    xL=xlim; yL=ylim;
    text(0.99*xL(2),0.99*yL(2), num2str(numTemporal_stationary),'HorizontalAlignment','right','VerticalAlignment','top')
    c = axis; line([0 0], [0 c(4)], 'LineStyle', '--', 'Color', 'r','LineWidth', 2)
    
    
    %% Saccade peth TOWARD PFD
    % 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 
    subtightplot(4,5,11,[Xtight Ytight]);
    
    
    
    
    %% Saccade peth AWAY FROM PFD
    % 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 
    subtightplot(4,5,11,[Xtight Ytight]);
    
    
    
    
    
    
    
    
    
    
    %% 
    subtightplot(4,5,20,[Xtight Ytight]);
    set(gca, 'YTick', [])
    title(sd.fn{iCell})
    
    if iCell < length(sd.fc)
        disp('press any key to continue')
        pause
    end
end

