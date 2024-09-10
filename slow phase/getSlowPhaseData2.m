function [data_out, data_outC, data_outM, data_outD, data_outMD] = getSlowPhaseData2(cfg_in)
% function [data_out, nSpikesRemoved] = getSlowPhaseData(cfg_in, sd, iCell)
%
%   Inputs
%           sd          - init file with sd.S
%           iCell       - which cell number to plot. leave empty [] if not calculating neural data
%           cfg_in      - config variable. Choose the "de-saccading" window.
%   Outputs
%           data_out    - the TSDs for this session, including horiz_eye_pos, horiz_eye_vel, horiz_eye_vel_smooth, AHV, and HD (a.k.a. "orientation")
%           data_outD   - the same TSDs, but with the pre and post saccade interval cut out (with antirestrict.m)
%           data_outMD  - same as above, but additionally with the stationary periods removed (with antirestrict.m)
%
%                           data_out.wheel                      - tsd of wheel data, in common timebase
%                           data_outM.wheel                     - tsd of wheel data, without stationary periods
%                           data_outD.wheel                     - tsd of wheel data, without saccade epochs
%                           data_outMD.wheel                    - tsd of wheel data, without stationary periods or saccade epochs
%                           data_outC.wheel                     - same as above, but with slow phase eye velocity cutoff 

%                           data_out.HD                         - tsd of HD (platform orientation), in common timebase
%                           data_outM.HD                        - tsd of HD (platform orientation), without stationary periods
%                           data_outD.HD                        - tsd of HD (platform orientation), without saccade epochs
%                           data_outMD.HD                       - tsd of HD (platform orientation), without stationary periods or saccade epochs
%                           data_outC.HD                        - same as above, but with slow phase eye velocity cutoff 

%                           data_out.AHV                        - tsd of AHV, in common timebase
%                           data_outM.AHV                       - tsd of AHV, without stationary epochs
%                           data_outD.AHV                       - tsd of AHV, without saccade epochs
%                           data_outMD.AHV                      - tsd of AHV, without stationary periods or saccade epochs
%                           data_outC.AHV                       - same as above, but with slow phase eye velocity cutoff 

%                           data_out.horiz_eye_pos              - tsd of EP, in common timebase
%                           data_outM.horiz_eye_pos             - tsd of EP, without stationary periods
%                           data_outD.horiz_eye_pos             - tsd of EP, without saccade epochs
%                           data_outMD.horiz_eye_pos            - tsd of EP, without stationary periods or saccade epochs
%                           data_outC.horiz_eye_pos             - same as above, but with slow phase eye velocity cutoff 

%                           data_out.horiz_eye_vel              - tsd of EV, in common timebase
%                           data_outM.horiz_eye_vel             - tsd of EV, without stationary periods
%                           data_outD.horiz_eye_vel             - tsd of EV, without saccade epochs
%                           data_outMD.horiz_eye_vel            - tsd of EV, without stationary periods or saccade epochs
%                           data_outC.horiz_eye_vel             - same as above, but with slow phase eye velocity cutoff 

%                           data_out.horiz_eye_vel_smooth       - tsd of filtered & smoothed EV, in common timebase
%                           data_outM.horiz_eye_vel_smooth      - tsd of filtered & smoothed EV, without stationary periods
%                           data_outD.horiz_eye_vel_smooth      - tsd of filtered & smoothed EV, without saccade epochs
%                           data_outMD.horiz_eye_vel_smooth     - tsd of filtered & smoothed EV, without stationary periods or saccade epochs
%                           data_outC.horiz_eye_vel_smooth      - same as above, but with slow phase eye velocity cutoff 


% Uses saccade .mat file to restrict data to slow phase only
% Useful for later tuning curve and GLM analysis
%
% 2024-04-18. MvdM, based on JJS getSlowPhaseTC3
% 2024-05-29. JJS. Added conditions to skip the resampling of neural data.
% 2024-09-04. JJS. Updated to use the variables that are saved in units of degrees and degrees/second (tsdHdeg, diffHdeg).
% 2024-09-08. JJS. This version (ver2) no longer takes sd or iCell as inputs. sd can take a long time to load, and there is no need to do anything with neurons at this step.
start = tic;
data_out = []; data_outD = []; data_outMD = [];

cfg_def = [];
cfg_def.FontSize = 20;
cfg_def.saccade_pre = 0.1;  % how many seconds to cut out before the saccade.
cfg_def.saccade_post = 0.1; % how many seconds to cut out after the saccade.
cfg_def.doPlotEYE = 0;
cfg_def.doPlotEVERYTHING = 0;
cfg_def.markerSize = 3;
cfg_def.plotWhichPhase = 'slow'; % 'slow', 'both'
cfg_def.medfilt_window_size = 11; % number of samples to use for slow phase velocity filter
cfg_def.smooth = 1;
cfg_def.doWarn = 1;
cfg_def.downsamplefactor = 10;
cfg_def.verbose = 1;
cfg_def.tightX = .04;
cfg_def.tightY = .04;
cfg_def.everythingPLOTfontsize = 15;
cfg_def.LineWidth = 1;
cfg_def.MarkerSize = 6;
cfg_def.restrict_velocity = 1; 
cfg_def.velocitythreshold = 35;  % velocity cutoff for eye velocity, in degrees per second. The restricted time periods are restricted from all data (HD,AHV,EP,EV)

cfg_Q = [];
cfg_Q.dt = 0.02; % binsize in s
cfg_Q.smooth = 'gauss'; % [], 'gauss'
cfg_Q.gausswin_size = 1; % gaussian window size in seconds; only used if cfg.smooth = 'gauss'
cfg_Q.gausswin_sd = 0.1; % SD for gaussian convolution in seconds
cfg_def.cfg_Q = cfg_Q;

cfg_master = ProcessConfig2(cfg_def, cfg_in);

%% Setup
SSN = HD_GetSSN; disp(SSN)
load(strcat(SSN, '-saccades-edited.mat'), 'AHV_tsd', 'combinedSaccades', 'tsdHdeg', 'diffHdeg') % load the edited (curated) saccade timestamps and the pupil trace time series (position and velocity)
load(strcat(SSN, '-AHV_StationaryTimes.mat'), 'STtstart', 'STtend')           % load the manually generated stationary period start/stop times

timestouse = ~isnan(combinedSaccades);                                      % create a logical of non-NaN saccade times. ***At some point, go through and remove the NaNs.
combinedSaccadesToUse = combinedSaccades(timestouse);                       % select only non-NaN saccade values to use

pre = combinedSaccadesToUse - cfg_master.saccade_pre;                       % get pre-saccade timestamps
post = combinedSaccadesToUse + cfg_master.saccade_post;                     % get post-saccade timestamps

if cfg_master.doWarn == 1; warning('off','all'); end                        % remove warnings, especially having to do with resampleTSD
tsdHdeg.data = tsdHdeg.data';                                               % .data needs to be 1 x n format to be a "well-formed tsd" for resampleTSD.m
diffHdeg.tvec = diffHdeg.tvec';                                             % .tvec needs to be n x 1 format to be a "well-formed tsd". *** Need to resave these vars so they have the correct (and consistent) format

%% RESAMPLE data to be on common timebase
cfg.method = 'rebinning';
new_tvec_edges = AHV_tsd.tvec(1): cfg_master.cfg_Q.dt: AHV_tsd.tvec(end); % take AHV_tsd from saccades-edited.mat file. This will start at t0 = 0 and tEnd = .tvec(end)
new_tvec = new_tvec_edges(1:end-1) + (cfg_master.cfg_Q.dt / 2);
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%       HEAD DIRECTION
tic;
if cfg_master.verbose == 1; disp('rebinning HD'); end
[~, orientation, starttimeUNIX, ~, ~, ~] = GetOrientationValues([]);
[orientationtouse] = downsampleOrientationValues(orientation, cfg_master.downsamplefactor); % if downsample factor = 10, it will go from 2kHz to 200Hz
head_direction_resamp = resampleTSD(cfg, orientationtouse, new_tvec);
toc;
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%       WHEEL
tic;
[wheel] = get_wheel_speed(starttimeUNIX);
if cfg_master.verbose == 1; disp('rebinning WHEEL'); end
wheel_resamp = resampleTSD(cfg, wheel, new_tvec);
toc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%        AHV
tic;
if cfg_master.verbose == 1; disp('rebinning AHV'); end
cfg = []; cfg.method = 'rebinning';
ahv_resamp = resampleTSD(cfg, AHV_tsd, new_tvec);
toc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%       EYE POSITION
tic;
if cfg_master.verbose == 1; disp('rebinning EP'); end
horiz_eye_pos_resamp = resampleTSD(cfg, tsdHdeg, new_tvec);
toc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%       EYE VELOCITY
tic;
if cfg_master.verbose == 1; disp('rebinning EV'); end
velocity_filt = medfilt1(diffHdeg.data, cfg_master.medfilt_window_size);                % apply median filter
velocity_smooth = smoothdata(velocity_filt);                                            % apply smoothing
horiz_eye_vel_smooth = tsd(diffHdeg.tvec, velocity_smooth);                             % make tsd
horiz_eye_vel_resamp = resampleTSD(cfg, diffHdeg, new_tvec);
horiz_eye_vel_smooth_resamp = resampleTSD(cfg, horiz_eye_vel_smooth, new_tvec);
toc

% Save resampled data
data_out.wheel = wheel_resamp;                                                          % WHEEL
data_out.HD = head_direction_resamp;                                                    % HD
data_out.AHV = ahv_resamp;                                                              % AHV
data_out.horiz_eye_pos = horiz_eye_pos_resamp;                                          % EP
data_out.horiz_eye_vel = horiz_eye_vel_resamp;                                          % EV
data_out.horiz_eye_vel_smooth = horiz_eye_vel_smooth_resamp;                            % EV, smoothed & filtered

%% Remove stationary periods only (M = "moving")
wheelM                  = antirestrict(wheel_resamp, STtstart, STtend);                 % WHEEL
head_directionM         = antirestrict(head_direction_resamp, STtstart, STtend);        % HD
ahvM                    = antirestrict(ahv_resamp, STtstart, STtend);                   % AHV
horiz_eye_posM          = antirestrict(horiz_eye_pos_resamp, STtstart, STtend);         % EP
horiz_eye_velM          = antirestrict(horiz_eye_vel_resamp, STtstart, STtend);         % EV
horiz_eye_vel_smoothM   = antirestrict(horiz_eye_vel_smooth_resamp, STtstart, STtend);  % EV, smoothed & filtered

% Save the resampled, MOVING-only data
data_outM.wheel                 = wheelM;                                               % WHEEL
data_outM.HD                    = head_directionM;                                      % HD
data_outM.AHV                   = ahvM;                                                 % AHV
data_outM.horiz_eye_pos         = horiz_eye_posM;                                       % EP
data_outM.horiz_eye_vel         = horiz_eye_velM;                                       % EV
data_outM.horiz_eye_vel_smooth  = horiz_eye_vel_smoothM;                                % EV, smooth & filtered

%% Remove quick phase periods only (D = "de-saccaded")
wheelD                  = antirestrict(wheel_resamp, pre, post);                        % WHEEL
head_directionD         = antirestrict(head_direction_resamp, pre, post);               % HD
ahvD                    = antirestrict(ahv_resamp, pre, post);                          % AHV
horiz_eye_posD          = antirestrict(horiz_eye_pos_resamp, pre, post);                % EP
horiz_eye_velD          = antirestrict(horiz_eye_vel_resamp, pre, post);                % EV
horiz_eye_vel_smoothD   = antirestrict(horiz_eye_vel_smooth_resamp, pre, post);         % EV, smoothed & filtered

% Save the resampled, DESACCADED data
data_outD.wheel                 = wheelD;                                               % WHEEL
data_outD.HD                    = head_directionD;                                      % HD
data_outD.AHV                   = ahvD;                                                 % AHV
data_outD.horiz_eye_pos         = horiz_eye_posD;                                       % EP
data_outD.horiz_eye_vel         = horiz_eye_velD;                                       % EV
data_outD.horiz_eye_vel_smooth  = horiz_eye_vel_smoothD;                                % EV, smooth & filtered

%% Remove quick phase periods from the MOVING-only data (M = "moving" + D = "de-saccaded")
wheelMD                  = antirestrict(wheelM, pre, post);                             % WHEEL
head_directionMD         = antirestrict(head_directionM, pre, post);                    % HD
ahvMD                    = antirestrict(ahvM, pre, post);                               % AHV
horiz_eye_posMD          = antirestrict(horiz_eye_posM, pre, post);                     % EP
horiz_eye_velMD          = antirestrict(horiz_eye_velM, pre, post);                     % EV
horiz_eye_vel_smoothMD   = antirestrict(horiz_eye_vel_smoothM, pre, post);              % EV, smoothed & filtered

% Save the resampled, DESACCADED data
data_outMD.wheel                 = wheelMD;                                             % WHEEL
data_outMD.HD                    = head_directionMD;                                    % HD
data_outMD.AHV                   = ahvMD;                                               % AHV
data_outMD.horiz_eye_pos         = horiz_eye_posMD;                                     % EP
data_outMD.horiz_eye_vel         = horiz_eye_velMD;                                     % EV
data_outMD.horiz_eye_vel_smooth  = horiz_eye_vel_smoothMD;                              % EV, smooth & filtered

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find velocity threshold times and restrict 
if cfg_master.restrict_velocity 
    i = abs(data_outMD.horiz_eye_vel_smooth.data) >= cfg_master.velocitythreshold;  % find the indices of moving, de-saccaded data that is ABOVE threshold (to remove later)
    
    W = data_outMD.wheel.data; W(i) = NaN; W = tsd(data_outMD.wheel.tvec, W);                                          % WHEEL
        
    HD = data_outMD.HD.data; HD(i) = NaN; HD = tsd(data_outMD.HD.tvec, HD);                                             % HD
    
    AHV = data_outMD.AHV.data; AHV(i) = NaN; AHV = tsd(data_outMD.AHV.tvec, AHV);                                       % AHV
    
    EP = data_outMD.horiz_eye_pos.data; EP(i) = NaN; EP = tsd(data_outMD.horiz_eye_pos.tvec, EP);                       % EP
    
    EV = data_outMD.horiz_eye_vel.data; EV(i) = NaN; EV = tsd(data_outMD.horiz_eye_vel.tvec, EV);                       % EV
    
    EVS = data_outMD.horiz_eye_vel_smooth.data; EVS(i) = NaN; EVS = tsd(data_outMD.horiz_eye_vel_smooth.tvec, EVS);     % EVS
end

% Save the resampled, MOVING-only data
data_outC.wheel                 = W;                % WHEEL
data_outC.HD                    = HD;               % HD
data_outC.AHV                   = AHV;              % AHV
data_outC.horiz_eye_pos         = EP;               % EP
data_outC.horiz_eye_vel         = EV;               % EV
data_outC.horiz_eye_vel_smooth  = EVS;              % EV, smooth & filtered

%% Plot EYE only
if cfg_master.doPlotEYE
    figure(1); clf;
    displayFactor = 10;
    plot(tsdHdeg.tvec, tsdHdeg.data .* displayFactor); hold on  % pupil position
    plot(diffHdeg.tvec, diffHdeg.data, 'r');  % pupil velocity
    plot(horiz_eye_vel_smooth.tvec, horiz_eye_vel_smooth.data, 'g', 'LineWidth', 2) % smoothed pupil velocity
    set(gca, 'FontSize', cfg_master.FontSize)
    plot(tsdHdeg.tvec, tsdHdeg.data .* displayFactor, '.b') % pupil position as points
    plot(diffHdeg.tvec, diffHdeg.data, 'r.') % pupil velocity as points
    plot(horiz_eye_vel_smooth.tvec, horiz_eye_vel_smooth.data, 'g', 'LineWidth', 5) % smoothed pupil velocity. Replotted here so the green line is in front.
    legend('pupil position', 'velocity', 'median filtered & smoothed velocity')
    title(SSN)
    xlabel('time (s)')
    ylabel('deg/s')
end

%% Plot EVERYTHING
if cfg_master.doPlotEVERYTHING
    figure(2); clf;
    %% WHEEL
    w1 = subtightplot(6,4,1, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_out.wheel.tvec, data_out.wheel.data, 'LineWidth', cfg_master.LineWidth)
    ylabel('WHEEL (cm/s)')
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    title('ALL DATA')
    set(gca, 'color', [0.9 0.9 0.9])
    
    w2 = subtightplot(6,4,2, [cfg_master.tightX cfg_master.tightY]);
    plot(data_outM.wheel.tvec, data_outM.wheel.data, 'r.', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    title('ROTATION ONLY')
    set(gca, 'color', [0.9 0.9 0.9])
    
    w3 = subtightplot(6,4,3, [cfg_master.tightX cfg_master.tightY]);
    plot(data_outD.wheel.tvec, data_outD.wheel.data, 'k.', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    title('DE-SACCADED')
    set(gca, 'color', [0.9 0.9 0.9])
    
    w4 = subtightplot(6,4,4, [cfg_master.tightX cfg_master.tightY]);
    plot(data_outMD.wheel.tvec, data_outMD.wheel.data, 'm.', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    title('ROTATION ONLY + DESACCADED')
    set(gca, 'color', [0.9 0.9 0.9])
    
    %% HD
    hd1 = subtightplot(6,4,5, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_out.HD.tvec, data_out.HD.data, 'LineWidth', cfg_master.LineWidth)
    ylabel('HD (deg)')
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    hd2 = subtightplot(6,4,6, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outM.HD.tvec, data_outM.HD.data, 'r.', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    hd3 = subtightplot(6,4,7, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outD.HD.tvec, data_outD.HD.data, 'k.', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    hd4 = subtightplot(6,4,8, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outMD.HD.tvec, data_outMD.HD.data, 'm.', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    %% AHV
    ahv1 = subtightplot(6,4,9, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_out.AHV.tvec, data_out.AHV.data, 'LineWidth', cfg_master.LineWidth)
    ylabel('AHV (deg/s)')
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    ahv2 = subtightplot(6,4,10, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outM.AHV.tvec, data_outM.AHV.data, 'r.', 'MarkerSize', cfg_master.MarkerSize);
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    ahv3 = subtightplot(6,4,11, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outD.AHV.tvec, data_outD.AHV.data, 'k.', 'MarkerSize', cfg_master.MarkerSize);
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    ahv4 = subtightplot(6,4,12, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outMD.AHV.tvec, data_outMD.AHV.data, 'm.', 'MarkerSize', cfg_master.MarkerSize);
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    %% EP
    ep1 = subtightplot(6,4,13, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_out.horiz_eye_pos.tvec, data_out.horiz_eye_pos.data, 'LineWidth', cfg_master.LineWidth)
    ylabel('EP (deg)')
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    ep2 = subtightplot(6,4,14, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outM.horiz_eye_pos.tvec, data_outM.horiz_eye_pos.data, 'r', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    ep3 = subtightplot(6,4,15, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outD.horiz_eye_pos.tvec, data_outD.horiz_eye_pos.data, 'k.', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    ep4 = subtightplot(6,4,16, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outMD.horiz_eye_pos.tvec, data_outMD.horiz_eye_pos.data, 'm.', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    %% EV
    ev1 = subtightplot(6,4,17, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_out.horiz_eye_vel.tvec, data_out.horiz_eye_vel.data, 'LineWidth', cfg_master.LineWidth)
    xlabel('time (s)')
    ylabel('EV (deg/s)')
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    ev2 = subtightplot(6,4,18, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outM.horiz_eye_vel.tvec, data_outM.horiz_eye_vel.data, 'r', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    ev3 = subtightplot(6,4,19, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outD.horiz_eye_vel.tvec, data_outD.horiz_eye_vel.data, 'k.', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    ev4 = subtightplot(6,4,20, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outMD.horiz_eye_vel.tvec, data_outMD.horiz_eye_vel.data, 'm.', 'MarkerSize', cfg_master.MarkerSize)
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    %% EV filtered & smoothed
    evf1 = subtightplot(6,4,21, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_out.horiz_eye_vel_smooth.tvec, data_out.horiz_eye_vel_smooth.data, 'LineWidth', cfg_master.LineWidth)
    xlabel('time (s)')
    ylabel('EV smoothed (deg/s)')
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])

    evf2 = subtightplot(6,4,22, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outM.horiz_eye_vel_smooth.tvec, data_outM.horiz_eye_vel_smooth.data, 'r.', 'MarkerSize', cfg_master.MarkerSize)
    xlabel('time (s)')
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    evf3 = subtightplot(6,4,23, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outD.horiz_eye_vel_smooth.tvec, data_outD.horiz_eye_vel_smooth.data, 'k.', 'MarkerSize', cfg_master.MarkerSize)
    xlabel('time (s)')
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    evf4 = subtightplot(6,4,24, [cfg_master.tightX cfg_master.tightY]); hold on
    plot(data_outMD.horiz_eye_vel_smooth.tvec, data_outMD.horiz_eye_vel_smooth.data, 'm.', 'MarkerSize', cfg_master.MarkerSize)
    xlabel('time (s)')
    set(gca, 'FontSize', cfg_master.everythingPLOTfontsize)
    set(gca, 'color', [0.9 0.9 0.9])
    
    linkaxes([w1 w2 w3 w4 hd1 hd2 hd3 hd4 ahv1 ahv2 ahv3 ahv4 ep1 ep2 ep3 ep4 ev1 ev2 ev3 ev4 evf1 evf2 evf3 evf4],'x')
    p1 = linkprop([w1 w2 w3 w4], 'ylim');
    p2 = linkprop([hd1 hd2 hd3 hd4], 'ylim');
    p3 = linkprop([ahv1 ahv2 ahv3 ahv4], 'ylim');
    p4 = linkprop([ep1 ep2 ep3 ep4], 'ylim');
    p5 = linkprop([ev1 ev2 ev3 ev4], 'ylim');
    p6 = linkprop([evf1 evf2 evf3 evf4], 'ylim');
end

if cfg_master.doWarn == 1; warning('on','all'); end  % Return warnings to ON state.
toc(start)

% clf
% plot(data_out.horiz_eye_vel.tvec, data_out.horiz_eye_vel.data); hold on; 
% axis([64 72 -200 200])
% plot(data_outMD.horiz_eye_vel_smooth.tvec, data_outMD.horiz_eye_vel_smooth.data, 'r.', 'MarkerSize', 25);
% plot(data_outC.horiz_eye_vel_smooth.tvec, data_outC.horiz_eye_vel_smooth.data, 'g.', 'MarkerSize', 15);
% set(gca, 'FontSize', 20)
% xlabel('time (s)')
% ylabel('eye vel (deg/s)')
% title(SSN)
% c = axis;
% line([c(1) c(2)], [35 35], 'Color', 'k', 'LineWidth', 3, 'LineStyle', '--')
% line([c(1) c(2)], [-35 -35], 'Color', 'k', 'LineWidth', 3, 'LineStyle', '--')
% legend('all data', 'desaccaded', 'desaccaded + slow phase cutoff')