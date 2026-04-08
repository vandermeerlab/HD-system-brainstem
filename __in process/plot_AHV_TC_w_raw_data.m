function plot_AHV_TC_w_raw_data(sd, iCell)

SSN = HD_GetSSN; disp(SSN);
FontSize = 18;
% Load Spikes
cfg = [];
cfg.uint = '64';
spikefiles = FindFiles('*.t');
cfg.fc = {spikefiles{iCell}};
Sold = LoadSpikes(cfg);

subtractStartTime = 1;
if subtractStartTime == 1 % New cheetah versions have timestamps
    if exist('events_ts.mat') == 2
        load('events_ts.mat');
    else
        events_ts = LoadEvents([]);
    end
    wrapper = @(events_ts) strcmp(events_ts, 'Starting Recording');
    A = cellfun(wrapper, events_ts.label);
    Startindex = find(A); % index which label says 'Start Recording'
    starttime = events_ts.t{Startindex}(1); % use the very first start record time
    for iC = 1:length(Sold.t)
        if strcmp(SSN, 'M054-2020-03-03') == 1
            S.t{iC} = Sold.t{iC} - Sold.t{iC}(1);
        else
            S.t{iC} = Sold.t{iC} - starttime;  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
        end
    end
end

myCell = SelectTS([], sd.S, iCell);
% get AHV Tuning Curve
cfg_AHV = [];
cfg_AHV.subsample_factor = 10;
[AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, myCell);
AHV_dt = median(diff(AHV_tsd.tvec));
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% #1 plot scatterplot
% p = subtightplot(6,6,1, [tightX tightY]);
cfg_Q = [];
cfg_Q.smooth = 'gauss';
cfg_Q.gausswin_sd = 0.05;
cfg_Q.dt = AHV_dt;
cfg_Q.tvec_edges = AHV_tsd.tvec(1):AHV_dt:AHV_tsd.tvec(end);
F = MakeQfromS(cfg_Q, myCell); % convert to FR
F.data = F.data ./ cfg_Q.dt;
F_idx = nearest_idx3(AHV_tsd.tvec, F.tvec);
AHV_F = F.data(:,F_idx);
ymax = max(AHV_F);

clf
set(gca, 'TickDir', 'out', 'FontSize', FontSize)
plot(AHV_tsd.data, AHV_F, '.', 'MarkerSize', .5, 'Color', [.5 .5 .5]); hold on
set(gca, 'Ylim', [0 ymax], 'FontSize', FontSize)
c = axis;
axis([-175 175 0 c(4)]);

LineWidth = 5;
plot(tc_out.binCenters, smoothdata(tc_out.tc), 'LineWidth', LineWidth, 'Color', 'k');




