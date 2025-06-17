function [out_n, out_t, peth_tvec, mean_FR_n, mean_FR_t, mn_shuff, mt_shuff, Z, w, percent_n, percent_t, nDir, tDir, cfg_out] = saccade_sig_distr_3(tfile, cfg_in)
% saccadePETHsig_ver1_1.m  This function determines whether each neuorn is significantly modulated around the saccade time
%
%   Inputs
%           tfile: a cell array in which each element is the file location for a given neuron to use. For example,
%                        {'C:\Jeff\U01\datatouse\M282\M282-2022-02-04-1\M282-2022-02-04-1-TT01_2.t'}
%                       These will be all of the neurons to analyze. Most likely, the list of Confidence = 1 NPH neurons.
%           cfg_in: config variable.

%   Outputs
%           a:
%           b:
%
%   To do:
%           (1) write a function to go through all of the saccade files and remove the NaNs and re-save
% ver1      This function calculates the distribution of maximum distances (i.e. differences in FR in a given bin) between
%           the true PETH and uniform distribution and shuffled PETH and uniform distribution.
% ver2      This version expands the subplot from a 2 x 2 to a 2 x 3 and adds plots of a larger peth window than what is used for signif. testing.
%           Also added option to save figure and save data to session folder.
% ver3      2024-09-16. This version works on a list of more than one neuron ... and arbitrarily long list of tfiles (as a cell array)

% warning('off', 'all')  % why isnt this working?

cfg_def.FontSize = 15;
cfg_def.numShuff = 1000; % how many shuffles to do
cfg_def.histbinnum = 500;
cfg_def.doPlot = 0;
cfg_def.peth_Window = [-.2 .1]; % the window for doing statistics on. These values should be smaller than cfg_def.window
cfg_def.window = [-1 1]; % the window for display
cfg_out = ProcessConfig(cfg_def, cfg_in);

if isempty(tfile)
    tfile = FindFiles('*.t');
end
numCells = length(tfile);
for iCell = 1:length(tfile)
    pushdir(fileparts(tfile{iCell}));
    temp = pwd; disp(temp)
    %% Load Spikes
    cfgS.uint = '64';
    cfgS.fc = {tfile{iCell}};
    myCell = LoadSpikes(cfgS);
    [csc_tsd, orientation, samplingrate, dt] = GetOrientationValues([]); %#ok<ASGLU>
    % fc = FindFiles('*.t', 'CheckSubdirs', 0);
    % [a, b, c] = fileparts(fc);
    % temp = strcmp(tfile, fc); % find which neuron we are talking about so it can be accessed later
    % neuronID = find(temp);
    
    %% Find and subtract Start Time
    % *********************** WRITE THIS AS ITS OWN FUNCTION *****************************
    dateswitch = datetime('2020-10-06');               % On this date I swtiched from using CSC21 to CSC33 for the platform encoder.
    SSN = HD_GetSSN;
    sessiondate = SSN(6:15);
    sessiondate = datetime(sessiondate);
    if sessiondate < dateswitch
        CSCtoUse = 21;     % CSC21 was used for the platform orientation encoder up until ______  % CSC33 was used after that date
    else
        CSCtoUse = 33;
    end
    cfg_csc = [];
    cfg_csc.fc = {FindFile(strcat('*CSC', num2str(CSCtoUse), '.ncs'))};
    cfg_csc.VoltageConvFactor = 10^6;
    csc_tsd = LoadCSC(cfg_csc);
    
    starttime = csc_tsd.tvec(1); % subtraction is necessary b/c cheetah stores timestamps in Linux time.
    endtime = csc_tsd.tvec(end);
    endtimetouse = endtime - starttime;
    myCell.t{1} = myCell.t{1} - starttime;  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
    aveFR = length(myCell.t{1})/endtimetouse;  % average firing rate for this neuron
    
    %% Load Saccade times
    load(FindFile('*saccades-edited.mat'), 'nasalSaccades', 'temporalSaccades', 'tsdH')
    keep = find(~isnan(temporalSaccades)); t = temporalSaccades(keep);  %#ok<*FNDSB>
    keep = find(~isnan(nasalSaccades)); n = nasalSaccades(keep);
    
    %% Calculate the true PETH
    cfg_MUA = [];
    %cfg_MUA.tvec = lfp_tsd.tvec'; % timebase to compute MUA on
    cfg_MUA.tvec = (0:0.001:endtimetouse); % timebase to compute MUA on
    cfg_MUA.tvec = cfg_MUA.tvec'; % flip it b/c tsd tvec needs to be n x 1, instead of 1 x n
    MUA_n = getMUA(cfg_MUA, myCell);
    MUA_t = getMUA(cfg_MUA, myCell);
    
    % warning('off', 'MATLAB:dispatcher:nameConflict')  % stop the annoying warnings from TSDpeth
    cfg_peth = [];
    cfg_peth.dt = 0.01;
    cfg_peth.doPlot = 0;
    cfg_peth.window = cfg_out.window; % this is the window for display purposes. the "peth_Window" is a smaller part of that for statistics.
    cfg_peth.mode = 'interp';
    out_n{iCell} = TSDpeth_fast(cfg_peth, MUA_n, n);
    out_t{iCell} = TSDpeth_fast(cfg_peth, MUA_t, t);
    
    peth_tvec(iCell,:) = out_n{iCell}.tvec; X.peth_tvec{iCell} = peth_tvec(iCell,:);
    statindicestouse = peth_tvec(iCell,:) >= cfg_out.peth_Window(1) & peth_tvec(iCell,:) <= cfg_out.peth_Window(2);   % what window do we want to use for doing statistics on? The full peth window is just for display.
    statwindow = peth_tvec(iCell,statindicestouse);
    X.statwindow = statwindow;
    nums = find(statindicestouse);
    firstbin = min(nums); windowstart = peth_tvec(iCell,firstbin);
    lastbin = max(nums); windowend = peth_tvec(iCell,lastbin);
    binnum = length(out_n{iCell}.tvec);
    binnumtouse = sum(statindicestouse); % sum up all the ones to get number of bins used for statistics
    stat_tvec = peth_tvec(iCell,statindicestouse);
    
    mean_FR_n(iCell) = mean(out_n{iCell}.data(statindicestouse));
    mean_FR_t(iCell) = mean(out_t{iCell}.data(statindicestouse));
    
    nMean_peth = repmat(mean_FR_n(iCell), 1, binnum);
    tMean_peth = repmat(mean_FR_t(iCell), 1, binnum);
    
    X.nDiff = abs(out_n{iCell}.data(statindicestouse) - nMean_peth(statindicestouse));
    X.tDiff = abs(out_t{iCell}.data(statindicestouse) - tMean_peth(statindicestouse));
    
    X.max_n = max(X.nDiff);
    X.max_t = max(X.tDiff);
    
    %% Caclulate the Circularly Shifted PETH
    tvec = tsdH.tvec; % tvec from the pupil position tsd
    mt_shuff{iCell} =  NaN(cfg_out.numShuff, binnumtouse);
    mn_shuff{iCell} = NaN(cfg_out.numShuff, binnumtouse);
    tic
    for iShuff = 1: cfg_out.numShuff
        h = waitbar(iShuff/cfg_out.numShuff);
        %     disp(iShuff);
        r(iShuff) = randsample(tvec,1); % choose a random value from the session times tvec
        nShift = n + r(iShuff); % shift nasal times by a random amount
        tShift = t + r(iShuff); % shift temporal times by a random amoun
        
        % SUBRACT SESSION DURATION FROM VALUES GREATER THAN THE LAST SESSION TIMESTAMP
        nList = nShift > endtimetouse; nIndices = find(nList);
        tList = tShift > endtimetouse; tIndices = find(tList); %
        
        nShift_new = nShift; tShift_new = tShift;
        if ~isempty(nIndices); nShift_new = nShift_new - endtimetouse; end  % subtract the session length from values > endtime
        if ~isempty(tIndices); tShift_new = tShift_new - endtimetouse; end  % subtract the session length from values > endtime
        
        % Calculate the shuffled PETH
        shuff_n = TSDpeth_fast(cfg_peth, MUA_n, nShift_new);
        mn_shuff{iCell}(iShuff,:) = shuff_n.data(statindicestouse);
        shuff_n_ave = nanmean(shuff_n.data(statindicestouse));
        shuff_n_uniform = repmat(shuff_n_ave, 1, binnum);
        
        shuff_t = TSDpeth_fast(cfg_peth, MUA_t, tShift_new);
        mt_shuff{iCell}(iShuff,:) = shuff_t.data(statindicestouse);
        shuff_t_ave = nanmean(shuff_t.data(statindicestouse));
        shuff_t_uniform = repmat(shuff_t_ave, 1, binnum);
        
        X.nDiff_shuff(iShuff,:) = abs(shuff_n.data(statindicestouse) - shuff_n_uniform(statindicestouse));  % *** not sure if I use true peth mean here or the shuffled peth mean.
        X.tDiff_shuff(iShuff,:) = abs(shuff_t.data(statindicestouse) - shuff_t_uniform(statindicestouse));
        
        X.max_n_shuff(iShuff) = max(X.nDiff_shuff(iShuff,:));
        X.max_t_shuff(iShuff) = max(X.tDiff_shuff(iShuff,:));
    end
    toc; disp('^^ time to calculate shuffled PETHs')
    
    set(groot,'ShowHiddenHandles','on')  % from matlab central
    delete(get(groot,'Children'))
    %% "Statistical test" - see where the sample value falls against the bootstrap distribution
    g_n = X.max_n_shuff < X.max_n;  % logical of how many elements in the distribution are less than the max distance
    h_n = nnz(g_n)/length(X.max_n_shuff);  % what is the fraction of distribution elements that are less than the max distance. e.g. 0.9 means the max distance is in the 90th percentile
    percent_n(iCell) = h_n;
    
    g_t = X.max_t_shuff < X.max_t;  % logical of how many elements in the distribution are less than the max distance
    h_t = nnz(g_t)/length(X.max_t_shuff);  % what is the fraction of distribution elements that are less than the max distance. e.g. 0.9 means the max distance is in the 90th percentile
    percent_t(iCell) = h_t;
    
    X.percent_n = percent_n(iCell);
    X.percent_t = percent_t(iCell);
    
    %% Determine whether the change in FR is an INCREASE or a DECREASE within the PETH window. 
    [max_n, max_index_n] = max(X.nDiff);
    bintouse = stat_tvec(max_index_n);
    [~, index_n] = min(abs(peth_tvec(iCell,:)-bintouse)); % find the bin in peth_tvec{iCell} which corresponds to the bin in stat_tvec that we know has the max difference
    nPeak = out_n{iCell}.data(index_n); 
    
    if nPeak > mean_FR_n(iCell) 
        nDir(iCell) = 1;  % Dir = direction of firing rate change. 1 is increased FR. -1 is decreased firing rate. 0 is equal. 
    elseif nPeak < mean_FR_n(iCell)
        nDir(iCell) = -1;
    elseif nPeak == mean_FR_n(iCell)
        nDir(iCell) = 0;
    else
        error('problem with determining firing rate change')
    end
    
    [max_t, max_index_t] = max(X.tDiff);
    bintouse = stat_tvec(max_index_t);
    [~, index_t] = min(abs(peth_tvec(iCell,:)-bintouse)); % find the bin in peth_tvec{iCell} which corresponds to the bin in stat_tvec that we know has the max difference
    tPeak = out_t{iCell}.data(index_t); 

    if tPeak > mean_FR_t(iCell) 
        tDir(iCell) = 1;  % Dir = direction of firing rate change. 1 is increased FR. -1 is decreased firing rate. 0 is equal. 
    elseif tPeak < mean_FR_t(iCell)
        tDir(iCell) = -1;
    elseif tPeak == mean_FR_t(iCell)
        tDir(iCell) = 0;
    else
        error('problem with determining firing rate change')
    end
    
    %% Plot It
    if cfg_out.doPlot == 1
        clf;
        subplot(321); % nasal peth
        plot(peth_tvec(iCell,:), out_n{iCell}.data);
        set(gca, 'FontSize', cfg_out.FontSize)
        ylabel('FR (Hz)')
        xlabel('time peri Saccade (s)'); c = axis;
        line([c(1) c(2)], [aveFR aveFR], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-'); c1 = axis;
        axis([cfg_out.window(1) cfg_out.window(2) 0 c1(4)])
        [~, cellID, ~] = fileparts(tfile{iCell});
        title(cellID)
        text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string', strcat('sess. ave FR = ', num2str(round(aveFR,1))))
        
        subplot(322) % temporal peth
        plot(peth_tvec(iCell,:), out_t{iCell}.data, 'Color', 'r');
        set(gca, 'FontSize', cfg_out.FontSize)
        xlabel('time peri Saccade (s)'); c = axis;
        line([c(1) c(2)], [aveFR aveFR], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-'); c2 = axis;
        axis([cfg_out.window(1) cfg_out.window(2) 0 c2(4)])
        title(cellID)
        text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string', strcat('sess. ave FR = ', num2str(round(aveFR,1))))
        z = max([c1(4) c2(4)]);
        axis([cfg_out.window(1) cfg_out.window(2) 0 z])
        line([windowstart windowstart], [0 z], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', [.5 .5 .5])   % start of statistics window
        line([windowend windowend], [0 z], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', [.5 .5 .5])   % end of statistics window
        legend('temporal', 'Location', 'Northwest')
        subplot(321);
        axis([cfg_out.window(1) cfg_out.window(2) 0 z])
        line([windowstart windowstart], [0 z], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', [.5 .5 .5])   % start of statistics window
        line([windowend windowend], [0 z], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', [.5 .5 .5])   % end of statistics window
        legend('nasal', 'Location', 'Northwest')
        
        subplot(323); % nasal peth
        plot(peth_tvec(iCell,:), out_n{iCell}.data);
        set(gca, 'FontSize', cfg_out.FontSize)
        ylabel('FR (Hz)')
        c3 = axis;
        axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) 0 c3(4)])
        line([c(1) c(2)], [mean_FR_n(iCell) mean_FR_n(iCell)], 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')   % aveFR line
        [a b] = max(X.nDiff);
        bintouse = stat_tvec(b);
        [val, index] = min(abs(peth_tvec(iCell,:)-bintouse)); % find the bin in peth_tvec{iCell} which corresponds to the bin in stat_tvec that we know has the max difference
        line([stat_tvec(b) stat_tvec(b)], [mean_FR_n(iCell) out_n{iCell}.data(index)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')   % plot a line showing the max diff btwn peth and aveFR
        
        
        subplot(324) % temporal peth
        plot(peth_tvec(iCell,:), out_t{iCell}.data, 'Color', 'r');
        set(gca, 'FontSize', cfg_out.FontSize)
        c4 = axis;
        axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) 0 c4(4)])
        line([c(1) c(2)], [mean_FR_t(iCell) mean_FR_t(iCell)], 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')   % aveFR line
        [a b] = max(X.tDiff);
        bintouse = stat_tvec(b);
        [val, index] = min(abs(peth_tvec(iCell,:)-bintouse)); % find the bin in peth_tvec{iCell} which corresponds to the bin in stat_tvec that we know has the max difference
        line([stat_tvec(b) stat_tvec(b)], [mean_FR_t(iCell) out_t{iCell}.data(index)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')   % plot a line showing the max diff btwn peth and aveFR
        axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) 0 z])
        subplot(323);
        axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) 0 z])
        
        subplot(325) %
        hist(X.max_n_shuff, cfg_out.histbinnum); c = axis;
        line([X.max_n  X.max_n], [c(3) c(4)], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '-')
        c5 = axis;
        xlabel('abs() max diff. of shuffles')
        ylabel('count')
        set(gca, 'FontSize', 15)
        title(strcat('percentile = ', num2str(percent_n*100)))
        text('Units', 'Normalized', 'Position', [0.5, 0.8], 'string', strcat('num shuff = ', num2str(cfg_out.numShuff)))
        
        subplot(326) %
        hist(X.max_t_shuff, cfg_out.histbinnum); c = axis;
        line([X.max_t  X.max_t], [c(3) c(4)], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '-')
        c6 = axis;
        xlabel('abs() max diff. of shuffles')
        ylabel('count')
        set(gca, 'FontSize', 15)
        z = max([c5(4) c6(4)]);
        axis([c6(1) c6(2) 0 z])
        line([X.max_t  X.max_t], [0 z], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '-')
        title(strcat('percentile = ', num2str(percent_t*100)))
        text('Units', 'Normalized', 'Position', [0.5, 0.8], 'string', strcat('num shuff = ', num2str(cfg_out.numShuff)))
        
        subplot(325)
        axis([c5(1) c5(2) 0 z])
        line([X.max_n  X.max_n], [0 z], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '-')
        disp('press any key to continue')
        pause
%         close
    end
    % warning('on', 'all')
    disp(strcat('percent nasal = ', num2str(X.percent_n)))
    disp(strcat('percent temporal = ', num2str(X.percent_t)))
    Z{iCell} = X;
    clear X
end

nUP = nDir == 1;                % cells with an increase in FR during nasal saccades, in the stats window 
nDOWN = nDir == -1;             % cells with an decrease in FR during nasal saccades, in the stats window 

tUP = nDir == 1;                % cells with an increase in FR during temporal saccades, in the stats window 
tDOWN = nDir == -1;             % cells with an increase in FR during temporal saccades, in the stats window 

w = [];
%% ONE OR THE OTHER significant, sign of change not considered
w.temporal_95 = percent_t >= 0.95; w.temporal_95_sum = sum(w.temporal_95); w.temporal_95_per = sum(w.temporal_95)/numCells;     % temporal sig. @ 95th percentile 
w.temporal_99 = percent_t >= 0.99; w.temporal_99_sum = sum(w.temporal_99); w.temporal_99_per = sum(w.temporal_99)/numCells;     % temporal sig. @ 99th percentile  

w.nasal_95 = percent_n >= 0.95; w.nasal_95_sum = sum(w.nasal_95); w.nasal_95_per = sum(w.nasal_95)/numCells;                    % nasal sig. @ 95th percentile 
w.nasal_99 = percent_n >= 0.99; w.nasal_99_sum = sum(w.nasal_99); w.nasal_99_per = sum(w.nasal_99)/numCells;                    % nasal sig. @ 95th percentile 

w.t_or_n_95 = percent_t >= 0.95 | percent_n >= 0.95; w.sum_t_or_n_95 = sum(w.t_or_n_95);  w.per_t_or_n_95 = w.sum_t_or_n_95/numCells;
w.t_or_n_99 = percent_t >= 0.99 | percent_n >= 0.99; w.sum_t_or_n_99 = sum(w.t_or_n_99);  w.per_t_or_n_99 = w.sum_t_or_n_99/numCells;

%% BOTH  significant (one up and one down) 
w.tup_and_ndown_95 = percent_t >= 0.95 & tDir == 1 &  percent_t >= 0.95 & nDir == -1;   % Increasing nasal & decreasing temporal @ 95th percentile
w.sum_tup_and_ndown_95 = sum(w.tup_and_ndown_95); w.per_tup_and_ndown_95 = w.sum_tup_and_ndown_95/numCells;

w.tdown_and_nup_95 = percent_t >= 0.95 & tDir == -1 &  percent_t >= 0.95 & nDir == 1;   % vice versa
w.sum_tdown_and_nup_95 = sum(w.tdown_and_nup_95); w.per_tdown_and_nup_95 = w.sum_tdown_and_nup_95/numCells;

w.tup_and_ndown_99 = percent_t >= 0.99 & tDir == 1 &  percent_t >= 0.99 & nDir == -1;   % Increasing nasal & decreasing temporal @ 99th percentile
w.sum_tup_and_ndown_99 = sum(w.tup_and_ndown_99); w.per_tup_and_ndown_99 = w.sum_tup_and_ndown_99/numCells;

w.tdown_and_nup_99 = percent_t >= 0.99 & tDir == -1 &  percent_t >= 0.99 & nDir == 1;   % vice versa 
w.sum_tdown_and_nup_99 = sum(w.tdown_and_nup_99); w.per_tdown_and_nup_99 = w.sum_tdown_and_nup_99/numCells; 

%% Significant Increasing Firing Rate Response  (either saccade type, or both) 
w.temporal_increase_95 = percent_t >= 0.95 & tDir == 1; w.sum_temporal_increase_95 = sum(w.temporal_increase_95); w.temporal_increase_95_percentage = w.sum_temporal_increase_95/numCells;
w.nasal_increase_95 = percent_n >= 0.95 & nDir == 1; w.sum_nasal_increase_95 = sum(w.nasal_increase_95); w.nasal_increase_95_percentage = w.sum_nasal_increase_95/numCells;

w.temporal_increase_99 = percent_t >= 0.99 & tDir == 1; w.sum_temporal_increase_99 = sum(w.temporal_increase_99); w.temporal_increase_99_percentage = w.sum_temporal_increase_99/numCells;
w.nasal_increase_99 = percent_n >= 0.99 & nDir == 1; w.sum_nasal_increase_99 = sum(w.nasal_increase_99); w.nasal_increase_99_percentage = w.sum_nasal_increase_99/numCells;

f95_t_incr = find(w.temporal_increase_95);  % temporal increase - find indices
f99_t_incr = find(w.temporal_increase_99); 

f95_n_incr = find(w.nasal_increase_95);  % nasal increase - find indices
f99_n_incr = find(w.nasal_increase_99); 

both_incr_95 = intersect(f95_t_incr, f95_n_incr); num_both_incr_95 = length(both_incr_95); per_both_95 = num_both_incr_95/numCells;
both_incr_99 = intersect(f99_t_incr, f99_n_incr); num_both_incr_99 = length(both_incr_99); per_both_99 = num_both_incr_99/numCells;

%% Significant Decreasing Firing Rate Response   (either saccade type, or both) 
w.temporal_decrease_95 = percent_t >= 0.95 & tDir == -1; w.sum_temporal_decrease_95 = sum(w.temporal_decrease_95); w.temporal_decrease_95_percentage = w.sum_temporal_decrease_95/numCells;
w.nasal_decrease_95 = percent_n >= 0.95 & nDir == -1; w.sum_nasal_decrease_95 = sum(w.nasal_decrease_95); w.nasal_decrease_95_percentage = w.sum_nasal_decrease_95/numCells;

w.temporal_decrease_99 = percent_t >= 0.99 & tDir == -1; w.sum_temporal_decrease_99 = sum(w.temporal_decrease_99); w.temporal_decrease_99_percentage = w.sum_temporal_decrease_99/numCells;
w.nasal_decrease_99 = percent_n >= 0.99 & nDir == -1; w.sum_nasal_decrease_99 = sum(w.nasal_decrease_99); w.nasal_decrease_99_percentage = w.sum_nasal_decrease_99/numCells;

w.f95_t_decr = find(w.temporal_decrease_95);  % temporal increase - find indices
w.f99_t_decr = find(w.temporal_decrease_99); 

w.f95_n_decr = find(w.nasal_decrease_95);  % nasal increase - find indices
w.f99_n_decr = find(w.nasal_decrease_99); 

w.both_decr_95 = intersect(w.f95_t_decr, w.f95_n_decr); w.num_both_decr_95 = length(w.both_decr_95); w.per_both_decr_95 = w.num_both_decr_95/numCells;
w.both_decr_99 = intersect(w.f99_t_decr, w.f99_n_decr); w.num_both_decr_99 = length(w.both_decr_99); w.per_both_decr_99 = w.num_both_decr_99/numCells;

%% NEITHER   (not significant for either saccade type)
w.numSig_neither_95 = percent_t < 0.95 & percent_n < 0.95; w.sumSig_neither_95 = sum(w.numSig_neither_95); w.per_Sig_neither_95 = w.sumSig_neither_95/numCells;
w.numSig_neither_99 = percent_t < 0.99 & percent_t < 0.99; w.sumSig_neither_99 = sum(w.numSig_neither_99); w.per_Sig_neither_99 = w.sumSig_neither_99/numCells;

w.neither_sig_95 = find(w.numSig_neither_95); w.num_neither_95 = length(w.neither_sig_95); w.per_neither_95 = w.num_neither_95/numCells;
w.neither_sig_99 = find(w.numSig_neither_99); w.num_neither_99 = length(w.neither_sig_99); w.per_neither_99 = w.num_neither_99/numCells;




% C = {'spn95', 'spt95', 'spn99', 'spt99', 'n_percent_95', 't_percent_95', 'n_percent_99', 't_percent_99', 'f95_n', 'f95_t', ...
%     'numSig_both', 'percent_sig_both', 'intersect', 'pn95', 'pt95', 'pn99', 'pt99'}; W = orderfields(w, C);

% desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
% Titles  = desktop.getClientTitles;
% for k = 1:numel(Titles)
%    Client = desktop.getClient(Titles(k));
%    if ~isempty(Client) & ...
%       strcmp(char(Client.getClass.getName), 'com.mathworks.mde.array.ArrayEditor')
%       Client.close();
%    end
% end