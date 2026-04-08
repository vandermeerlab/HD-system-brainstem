function [] = saccadeTuning(varargin)

% 4/2021. JJS.
% For each neuron, determine whether FR is signif. different from shuffled spikes

process_varargin(varargin);

fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    cfg_in.doPlotThresholds = 0;
    [temporalSaccades, nasalSaccades, combinedSaccades, index_tP_final, index_nP_final, tsdH, tsdV, diffH, diffV] = processPupilData2(cfg_in); %#ok<ASGLU> % get saccade times
    
    [S] = LoadSpikesJeff;
    for iCell = 1:length(S.t)
        
    t = temporalSaccades;
    cfg_in = [];
    cfg_in.window = [-.5 .5];
    cfg_in.dt = 0.005;
    [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm(cfg_in, S, t); 
    
    m = histc(outputS, outputIT);
    FRxBin = m/cfg.dt/length(t);
%     bar(outputIT,FRxBin);
    d = diff(S.t{1});               % get the ISIs
    r = randperm(length(d));        % permute the order
    R = ts(r);
    [outputS, outputT, outputGau, outputIT, cfg] = SpikePETHvdm(cfg_in, S, t); 
    
%     cumsum(randperm(diff(S.t{1})))         % why cumsum?
    
    end
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    s_idx = find(S.t{iCell} > temporalSaccades -0.5 & S.t{iCell} < temporalSaccades +0.5);