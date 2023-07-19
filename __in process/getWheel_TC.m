function [tc_wheel] = getWheel_TC(cfg_in, sd)
% JJS. 2023-05-17.
% Calculate and plot the tuning curve for wheel speed
% input:   sd - session data structure with spike trains S
% output:  pupilTC - a 1 x nCell array of structures, each with the fields tc, occ_hist, spk_hist, usr, cfg

doPlot = 1;
cfg_def.FontSize = 16;
cfg_def.speedthresh = 0.3; % cm / sec  [minimum wheel speed for TC]
cfg_def.occthresh = 0.5;  % seconds [minimum occupancy for each speed bin]
cfg = ProcessConfig2(cfg_def, cfg_in);


updownTSD = getQEupdown([]);
state_tsd = ConvertQEUpDownToState(updownTSD);
[~, wheel_tsd] = ConvertQEStatesToAngle([], state_tsd);
[d, speed, cfg_w] = ConvertWheeltoSpeed([], wheel_tsd);   % d = distance.  speed = speed of the wheel in cm/sec
speed.tvec = speed.tvec - sd.starttime;
speed.data = -speed.data; % we want forward motion to be displayed as a positive velocity

% calculate Q matrix
speed_dt = median(diff(speed.tvec));  % tsdH is the pupil position variable
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = speed_dt;
cfg_Q.tvec_edges = speed.tvec(1): speed_dt: speed.tvec(end);

%
F = MakeQfromS(cfg_Q, sd.S); % convert to FR
F.data = F.data ./ cfg_Q.dt;
F_idx = nearest_idx3(speed.tvec, F.tvec); % find FR corresponding to each AHV sample
tsd_W = F.data(:,F_idx);

cfg_tc = [];
cfg_tc.nBins = 50;
cfg_tc.binEdges = {linspace(-5, 30, 101)};
cfg_tc.occ_dt = median(diff(speed.tvec));
cfg_tc.minOcc = 1;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
tc_wheel = TuningCurves(cfg_tc, sd.S, speed);

if doPlot == 1
    for iCell = 1:length(sd.S.t)
        clf
        yyaxis right
        z = speed.data > cfg.speedthresh | speed.data < -cfg.speedthresh;
        plot(speed.data(z), tsd_W(iCell,z), '.', 'MarkerSize', .5, 'color', [.8 .8 .8]); hold on
        %         h = lsline;
        %         set(h(1), 'Color', 'k')
        %         set(h(1), 'LineWidth', 2)
        %         speeddata = speed.data(z)';
        %         speeddata(:,2) = ones(length(speeddata), 1);
        %         [~,~,~,~,stats] = regress(tsdH_F{iCell}(1,z)', speeddata);
        c = axis;
        line([-cfg.speedthresh -cfg.speedthresh], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
        line([cfg.speedthresh cfg.speedthresh], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
        %         text(NaN, NaN, strcat('Rsq =', sprintf('%0.2f', stats(1))), 'FontSize', 12, 'Units', 'normalized', 'Position', [.55 .85 0])
        
        plot(tc_wheel.binCenters(tc_wheel.occ_hist > cfg.occthresh), smoothdata(tc_wheel.tc(iCell,(tc_wheel.occ_hist > cfg.occthresh))), 'k', 'LineWidth', 3);
        axis tight
        ylabel('FR (Hz)', 'FontSize', cfg.FontSize)
        xlabel('Speed (cm/s)', 'FontSize', cfg.FontSize)
        title('Wheel Speed Tuning Curve', 'FontSize', cfg.FontSize)
        yyaxis left
        set(gca, 'YTick', [])
        yyaxis right
        grid on
        set(gca, 'FontSize', cfg.FontSize)
        p.XAxisLocation = 'top';
        axis tight
        if iCell ~= length(sd.S.t)
            disp('press any key to continue')
            pause
        end
    end
end