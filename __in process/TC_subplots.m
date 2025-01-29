function TC_subplots(sd)

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
    subtightplot(3,3,1,[Xtight Ytight]); hold on
    %     set(gca, 'YTick', [])
    set(gca, 'XTick', -180:30:180)
    %     grid on
    ax = gca;
    ax.XGrid = 'on';
    ax.GridLineStyle = '-';
    
    title('HD TC', 'FontSize', FontSizeToUse)
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
    
    
    %% AHV TC
    subtightplot(3,3,2,[Xtight Ytight]);
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
    title('AHV TC', 'FontSize', FontSizeToUse)
    
    %% Wheel TC
    subtightplot(3,3,3,[Xtight Ytight]);
    set(gca, 'YTick', [])
    
    tc_wheel = TuningCurves(cfg_wheel, myCell, speed);
    plot(tc_wheel.binCenters, tc_wheel.tc, 'k', 'LineWidth', 3); hold on
    plot(tc_wheel.binCenters, smoothdata(tc_wheel.tc), 'LineWidth', 3, 'Color', 'r', 'LineStyle', '--');
    
    ax = gca;
    ax.XGrid = 'on';
    ax.GridLineStyle = '-';
    title('Wheel TC', 'FontSize', FontSizeToUse)
    
    %% HD Occupancy
    subtightplot(3,3,4,[Xtight Ytight]);
    plot(tc_HD.binCenters, tc_HD.occ_hist)
    %     set(gca, 'YTick', [])
    set(gca, 'XTick', -180:30:180)
    grid on
    title('HD Occupancy', 'FontSize', FontSizeToUse)
    ylabel('Seconds', 'FontSize', FontSizeToUse)
    
    %% AHV Occupancy
    subtightplot(3,3,5,[Xtight Ytight]);
    plot(tc_AHV.binCenters, tc_AHV.occ_hist)
    c = axis;
    axis([c(1) c(2) 0 15]);
    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    set(gca, 'YTick', 0:3:15)
    title('AHV Occupancy', 'FontSize', FontSizeToUse)
    
    %% Wheel Occupancy
    subtightplot(3,3,6,[Xtight Ytight]);
    plot(tc_wheel.binCenters, tc_wheel.occ_hist)
    c = axis;
    axis([c(1) c(2) 0 15]);
    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    ax.GridLineStyle = '-';
    set(gca, 'YTick', 0:3:15)
    title('Wheel Occupancy', 'FontSize', FontSizeToUse)
    
    %% Pupil Position TC 
    subtightplot(3,3,7,[Xtight Ytight]);
    set(gca, 'YTick', [])
    
    cfg_pupil.binEdges = {linspace(-60, 60, 61)};
    cfg_pupil.occ_dt = median(diff(sd.tsdH.tvec));
    cfg_pupil.minOcc = 10;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
    tc_pupil_position = TuningCurves([], myCell, sd.tsdH);
    
    plot(tc_pupil_position.binCenters, tc_pupil_position.tc, 'LineWidth', 3, 'Color', 'k'); hold on
    plot(tc_pupil_position.binCenters, smoothdata(tc_pupil_position.tc), 'LineWidth', 3, 'Color', 'r', 'LineStyle', '--');
    c = axis;
    axis([-30 30 0 c(4)]);
    % axis tight
    
    title('Pupil Position TC', 'FontSize', FontSizeToUse)
    ylabel('FR (Hz)', 'FontSize', FontSizeToUse)
    xlabel('Pixels')
    
    %% Pupil Velocity TC 
    subtightplot(3,3,8,[Xtight Ytight]);
    set(gca, 'YTick', [])
    cfg_pupil_velocity = [];
    cfg_pupil_velocity.binEdges = {linspace(-40, 40, 41)};
    tc_pupil_velocity.occ_dt = median(diff(sd.diffH.tvec));
    tc_pupil_velocity.minOcc = 100;  % remember that Occ is measured in samples, not in seconds. Usually 20ms per sample, b/c the camera sampling rate is 50Hz,
    
    tc_pupil_velocity = TuningCurves(cfg_pupil_velocity, myCell, sd.diffH);
    plot(tc_pupil_velocity.binCenters, tc_pupil_velocity.tc, 'LineWidth', 3, 'Color', 'k'); hold on
    plot(tc_pupil_velocity.binCenters, smoothdata(tc_pupil_velocity.tc), 'LineWidth', 3, 'Color', 'r', 'LineStyle', '--');
    xlabel('Pixels/sec')
    title('Pupil Velocity TC', 'FontSize', FontSizeToUse)
    
    c = axis;
    axis([-30 30 0 c(4)]);
    %% 
    subtightplot(3,3,9,[Xtight Ytight]);
    set(gca, 'YTick', [])
    title(sd.fn{iCell})
    
    if iCell < length(sd.fc)
        disp('press any key to continue')
        pause
    end
end

