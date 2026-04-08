function [] = plotStationaryVsSpontaneous

% JJS. 2022-10-18.
% This function plots peths for saccades when the platform was moving (i.e. 'evoked saccades') and peths for saccades when the platform was stationary ('spontaneous').
% 2 x 4 subplot of peths, tuning curves, and firing rate autocorrelation

occthresh = 0.5;
LineWidth = 4;
FontSize = 12;
% cfg_in = [];
[~, ~, ~, ~, nasal_timestamps_MOVING, temporal_timestamps_MOVING] = isolateManualSaccades();
[~, ~, ~, ~, nasal_timestamps_REST, temporal_timestamps_REST] = isolateStationarySaccades();
close all;

S = LoadSpikesJeff;

for iCell = 1:length(S.t)
    disp(strcat('myCell', num2str(iCell)))
    myCell = SelectTS([], S, iCell);
    figure(iCell)
    
    numPanels = 8;
    for iPanel = 1:numPanels
        switch iPanel
            case 1
                %% plot nasal and temporal saccasde peths MOVING, -200 ms to +200ms
                subplot(2,4,1);
                %                 cfg_1.doRaster = 0;
                %                 cfg_1.doBar = 1;
                cfg_1.doPlot = 0;
                cfg_1.window = [-2 2];
                cfg_1.dt = 0.05;
                [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, myCell, nasal_timestamps_MOVING);
                [mn1, edges] = histcounts(outputS_n, outputIT_n);
                plot(edges(1:end-1), mn1/cfg_1.dt/length(nasal_timestamps_MOVING), 'LineWidth', LineWidth); % this is a hack. should replace binedgges with bincenters
                
                hold on
                [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, myCell, temporal_timestamps_MOVING);
                [mt1, edges] = histcounts(outputS_t, outputIT_t);
                %                 Z = mt1/cfg_1.dt/length(temporal_timestamps_MOVING);
                plot(edges(1:end-1), mt1/cfg_1.dt/length(temporal_timestamps_MOVING), 'LineWidth', LineWidth);
                title('Saccade peth MOVING')
                c = axis;
                line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
                legend('nasal', 'temporal', '', 'FontSize', FontSize)
                set(gca, 'FontSize', FontSize)
                ylabel('Firing Rate (Hz)')
                xlabel('time peri saccade (s)')
            case 2
                %% plot nasal and temporal saccasde peths STATIONARY, -200 ms to +200ms
                subplot(2,4,2);
                cfg_1.doPlot = 0;
                cfg_1.window = [-2 2];
                cfg_1.dt = 0.005;
                [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, myCell, nasal_timestamps_REST);
                [mn2, edges] = histcounts(outputS_n, outputIT_n);
                plot(edges(1:end-1), mn2/cfg_1.dt/length(nasal_timestamps_MOVING), 'LineWidth', LineWidth);
                
                hold on
                [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, myCell, temporal_timestamps_REST);
                [mt2, edges] = histcounts(outputS_t, outputIT_t);
                plot(edges(1:end-1), smoothdata(mt2/cfg_1.dt/length(temporal_timestamps_MOVING)), 'LineWidth', LineWidth);
                title('Saccade peth STATIONARY')
                c = axis;
                line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
                legend('nasal', 'temporal', '', 'FontSize', FontSize)
                set(gca, 'FontSize', FontSize)
                ylabel('Firing Rate (Hz)')
                xlabel('time peri saccade (s)')
            case 3
                %% plot nasal and temporal saccasde peths MOVING, -2 seconds to +2 seconds
                subplot(2,4,5);
                cfg_1.doPlot = 0;
                cfg_1.window = [-.2 .2];
                cfg_1.dt = 0.01;
                [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, myCell, nasal_timestamps_MOVING);
                [mn3, edges] = histcounts(outputS_n, outputIT_n);
                plot(edges(1:end-1), mn3/cfg_1.dt/length(nasal_timestamps_MOVING), 'LineWidth', LineWidth);
                
                hold on
                [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, myCell, temporal_timestamps_MOVING);
                [mt3, edges] = histcounts(outputS_t, outputIT_t);
                plot(edges(1:end-1), mt3/cfg_1.dt/length(temporal_timestamps_MOVING), 'LineWidth', LineWidth);
                title('Saccade peth MOVING: zoomed in')
                c = axis;
                line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
                legend('nasal', 'temporal', '', 'FontSize', FontSize)
                ylabel('Firing Rate (Hz)')
                xlabel('time peri saccade (s)')
            case 4
                %% plot nasal and temporal saccasde peths MOVING, -2 seconds to +2 seconds
                subplot(2,4,6);
                cfg_1.doPlot = 0;
                cfg_1.window = [-.2 .2];
                cfg_1.dt = 0.01;
                [outputS_n, ~, ~, outputIT_n, ~] = SpikePETH_either(cfg_1, myCell, nasal_timestamps_REST);
                [mn4, edges] = histcounts(outputS_n, outputIT_n);
                plot(edges(1:end-1), mn4/cfg_1.dt/length(nasal_timestamps_MOVING), 'LineWidth', LineWidth);
                
                hold on
                [outputS_t, ~, ~, outputIT_t, ~] = SpikePETH_either(cfg_1, myCell, temporal_timestamps_REST);
                [mt4, edges] = histcounts(outputS_t, outputIT_t);
                plot(edges(1:end-1), mt4/cfg_1.dt/length(temporal_timestamps_MOVING), 'LineWidth', LineWidth);
                title('Saccade peth STATIONARY: zoomed in')
                c = axis;
                line([0 0], [c(3) c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
                legend('nasal', 'temporal', '', 'FontSize', FontSize)
                set(gca, 'FontSize', FontSize)
                ylabel('Firing Rate (Hz)')
                xlabel('time peri saccade (s)')
            case 5
                
                
            case 6
                
                
                
            case 7
                %% Pupil position tuning curve
                subplot(2,4,7)
                try 
                    load(FindFile('*saccades-edited.mat'), 'tsdH')  
                catch
                    warning('cannot find saccade data')
                end
                % tsdH is the pupil position variable
                tsdH.data = tsdH.data'; % change the shape so that it is a "well-formed tsd" for tuning curves
                cfg_tc = [];
                cfg_tc.nBins = 50;
                cfg_tc.binEdges = {linspace(-60, 60, 101)};
                cfg_tc.occ_dt = median(diff(tsdH.tvec));
                cfg_tc.minOcc = 10;  % remember that Occ is measured in samples (usually 5ms per sample), not in seconds
                tc_pupil = TuningCurves(cfg_tc, S, tsdH);
                plot(tc_pupil.usr.binCenters(tc_pupil.occ_hist>occthresh), smoothdata(tc_pupil.tc(iCell,(tc_pupil.occ_hist>occthresh))), 'k', 'LineWidth', 3);
                set(gca, 'FontSize', FontSize)
                ylabel('Firing Rate (Hz)')
                xlabel('Pupil Position (pixels)')
                
            case 8
                %% AHV tuning curve
                subplot(2,4,8)
                % calculate AHV and get tuning curve
                cfg_AHV = [];
                cfg_AHV.subsample_factor = 10;
                [AHV_tsd, tc_AHV] = AHV_tuning(cfg_AHV, S);
                AHV_dt = median(diff(AHV_tsd.tvec));
                % calculate Q matrix
                cfg_Q = [];
                cfg_Q.smooth = 'gauss';
                cfg_Q.gausswin_sd = 0.05;
                cfg_Q.dt = AHV_dt;
                cfg_Q.tvec_edges = AHV_tsd.tvec(1):AHV_dt:AHV_tsd.tvec(end);
                F = MakeQfromS(cfg_Q, S); % convert to FR
                % convert to FR
                F.data = F.data ./ cfg_Q.dt;
                % find FR corresponding to each AHV sample
                F_idx = nearest_idx3(AHV_tsd.tvec, F.tvec);
                AHV_F = F.data(:,F_idx);
                plot(AHV_tsd.data, AHV_F(iCell,:), '.', 'MarkerSize', .5, 'color', [.8 .8 .8]);
                hold on
                plot(tc_AHV.usr.binCenters(tc_AHV.occ_hist>occthresh), smoothdata(tc_AHV.tc(iCell,(tc_AHV.occ_hist>occthresh))), 'k', 'LineWidth', 3);
                set(gca, 'FontSize', FontSize)
                xlabel('AHV (deg/s)')
                ylabel('Firing Rate (Hz)')
        end
    end
    disp('press any key to continue')
    pause
end



