function FindSaccades(varargin)

doPlot = 1;
% CheckFilteredSignal = 1;
FontSize = 20;
c =  [355  390  -25   25];
threshH = 15;
threshL = -15;
process_varargin(varargin)

% SSN = HD_GetSSN;
f = FindFiles('*VT1_proc.mat');
load(f{1}, 'pupil');

%% Get TimeStamps
[~, videofn, ext] = fileparts(FindFiles('*VT1.nvt'));  % get timestamps from the video tracking data
cfg = [];
cfg.fn = strcat(videofn, ext);
cfg.removeZeros = 0 ;
pos_tsd = LoadPos(cfg);
pupiltime = pos_tsd.tvec;   % it apprears to Nvt file is 2 frames longer than the number of frames from facemap. ask matt about this
pupiltime = pupiltime - pupiltime(1);


%% Create TSD
pupildata = pupil{1}.com(:,2) - mean(pupil{1}.com(:,2));
pupiltsd = tsd(pupiltime(2:end-1),pupildata');

%% Filter the signal
cfg_filter = [];
cfg_filter.type = 'butter'; 
cfg_filter.band = '';
cfg_filter.f = [2 500];
pupiltsd.cfg.hdr{1}.Fs = 1 ./ median(diff(pupiltsd.tvec));   % append the sampling rate
pupil_tsdF = FilterLFP(cfg_filter, pupiltsd);

diffp = diff(horiz_tsd.data);
diffp_tsd = tsd(horiz_tsd.tvec(2:end)', diffp); 
diffp_tsd.cfg.hdr{1}.Fs = 1 ./ median(diff(horiz_tsd.tvec)); 
diffp_tsdF = FilterLFP(cfg_filter, diffp_tsd);


H = find(diffp_tsd.data > threshH);
L = find(diffp_tsd.data < threshL);












tstart = horiz_tsd.tvec(1);
tend = horiz_tsd.tvec(end); 
if doPlot == 1
    clf
    a1 = subplot(3,1,1);
    plot(horiz_tsd.tvec, horiz_tsd.data);         % pupil{1}.com(:,2)  is the horizontal eye position from facemap
    ylabel('Horiz. Camera axis (pixels)', 'FontSize', FontSize)
    xlabel('Time (sec)', 'FontSize', FontSize)
    set(gca, 'FontSize', FontSize)
    axis([c(1) c(2) c(3) c(4)])
    
    a2 = subplot(3,1,2);
    plot(horiz_tsdF.tvec, horiz_tsdF.data)
    ylabel('diff  (pixels)', 'FontSize', FontSize)
    xlabel('Time (sec)', 'FontSize', FontSize)
    set(gca, 'FontSize', FontSize)
    axis([c(1) c(2) c(3) c(4)])
    
    a3 = subplot(3,1,3); 
    hold on
    plot(diffp_tsd.tvec, diffp_tsd.data)
    ylabel('diff  (pixels)', 'FontSize', FontSize)
    xlabel('Time (sec)', 'FontSize', FontSize)
    set(gca, 'FontSize', FontSize)
    axis([c(1) c(2) c(3) c(4)])
    
    line([tstart tend], [threshL threshL], 'Color', 'k')
    line([tstart tend], [threshH threshH], 'Color', 'k')
    plot(diffp_tsd.tvec(H), diffp_tsd.data(H), 'm.', 'MarkerSize', 10)
    plot(diffp_tsd.tvec(L), diffp_tsd.data(L), 'g.', 'MarkerSize', 10)
    
    
    linkaxes([a1 a2 a3])
%     subplot(numPlots,1,4)
%     plot(diffp_tsdF.tvec, yy)
%     ylabel('diff  (pixels)', 'FontSize', FontSize)
%     xlabel('Time (sec)', 'FontSize', FontSize)
%     axis([c(1) c(2) c(3) c(4)])
end















