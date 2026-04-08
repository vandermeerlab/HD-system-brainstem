function [TC_all, TC_norm, TC_cellID] = get_TuningCurves(cfg_in, fd)

subtractStartTime = 1;
cfg_def = [];
cfg = ProcessConfig(cfg_def,cfg_in);

TC_all.tc = [];
TC_all.bins = [];
TC_cellID = [];

if isempty(fd)
    fd = FindFiles('*keys.m');
end
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN;
    disp(SSN)
    
    cfg_S = [];
    cfg_S.uint = '64';
    S = LoadSpikes(cfg_S);
    
    if subtractStartTime == 1 % New cheetah versions have timestamps
        events_ts = LoadEvents([]);
        assert(strcmp(events_ts.label{1}, 'Starting Recording'))=1;
        for iC = 1:length(S.t)
            S.t{iC} = S.t{iC} - events_ts.t{1}(1);  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
        end
    end
    % default number of bins here is 100
    cfg_AHV = [];
    cfg_AHV.subsample_factor = 10;
    [AHV_tsd, tc_out] = AHV_tuning(cfg_AHV, S);
    numCells = size(tc_out.tc, 1);
    %     plot(tc_out.usr.binCenters, tc_out.tc(iCell,:), 'LineWidth', LineWidth);
    TC_all.tc = vertcat(TC_all.tc, tc_out.tc);
    TC_all.bins = vertcat(TC_all.bins, repmat(tc_out.binCenters, numCells, 1));
    clear tc_out
    
    TC_cellID = vertcat(TC_cellID, S.label');
    popdir;
end

TC_norm.tc = [];
TC_norm.bins = TC_all.bins;
for iCell = 1:size(TC_all.tc,1)
    TC_norm.tc(iCell,:) = TC_all.tc(iCell,:)./max(TC_all.tc(iCell,:));
end

% [B, I] = max(TC_norm.tc');


