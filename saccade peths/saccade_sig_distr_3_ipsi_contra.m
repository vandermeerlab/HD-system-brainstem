function [out_ipsi, out_contra, peth_tvec, mean_FR_ipsi, mean_FR_contra, ipsi_shuff, contra_shuff, Z, w, percent_ipsi, percent_contra, ipsi_Dir, contra_Dir, cfg_out] = saccade_sig_distr_3_ipsi_contra(tfile, cfg_in)
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
%               (1) write a function to go through all of the saccade files and remove the NaNs and re-save
% ver1          This function calculates the distribution of maximum distances (i.e. differences in FR in a given bin) between
%               the true PETH and uniform distribution and shuffled PETH and uniform distribution.
% ver2          This version expands the subplot from a 2 x 2 to a 2 x 3 and adds plots of a larger peth window than what is used for signif. testing.
%               Also added option to save figure and save data to session folder.
% ver3          2024-09-16. This version works on a list of more than one neuron ... and arbitrarily long list of tfiles (as a cell array)

% ipsi-contra   2024-09-20. This version uses ipsi-contra saccade variables instead of temporal and nasal saccades.  

%   Key to variable naming. ALWAYS FOR THE LEFT EYE! 

%    LEFT hemisphere neuron, temporal (leftward) saccade   = IPSI 
%    LEFT hemisphere neuron, nasal (rightward) saccade      = CONTRA 

%    RIGHT hemisphere neuron, temporal (leftward) saccade  = CONTRA 
%    RIGHT hemisphere neuron, nasal (rightward) saccade     = IPSI 

% warning('off', 'all')  % why isnt this working? for the CheckTSD warnings 

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
    EvalKeys;                           % load keys 
    
    %% Load Spikes
    cfgS.uint = '64';
    cfgS.fc = {tfile{iCell}};
    myCell = LoadSpikes(cfgS);
    [csc_tsd, orientation, samplingrate, dt] = GetOrientationValues([]); %#ok<ASGLU>

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
    
    %% Assign IPSI and CONTRA saccades for this individual neuron 
    tdSession = FindFiles('*.t'); % how many neurons in this session
    neuronArray = strcmp(tfile{iCell}, tdSession);
    iNeuron = find(neuronArray); 
    if strcmp(ExpKeys.Hemisphere{iNeuron}, 'L')  
        ipsi = t; contra = n; 
    elseif strcmp(ExpKeys.Hemisphere{iNeuron}, 'R')
        ipsi = n; contra = t; 
    else
        error('this neuron does not have a correct hemsiphere designation')   % what about center {'C'} cells? I think I removed them. 
    end
    
    %% Calculate the true PETH
    cfg_MUA = [];
    %cfg_MUA.tvec = lfp_tsd.tvec'; % timebase to compute MUA on
    cfg_MUA.tvec = (0:0.001:endtimetouse); % timebase to compute MUA on
    cfg_MUA.tvec = cfg_MUA.tvec'; % flip it b/c tsd tvec needs to be n x 1, instead of 1 x n
    MUA = getMUA(cfg_MUA, myCell);
    
    % warning('off', 'MATLAB:dispatcher:nameConflict')  % stop the annoying warnings from TSDpeth
    cfg_peth = [];
    cfg_peth.dt = 0.01;
    cfg_peth.doPlot = 0;
    cfg_peth.window = cfg_out.window; % this is the window for display purposes. the "peth_Window" is a smaller part of that for statistics.
    cfg_peth.mode = 'interp';
    out_ipsi{iCell} = TSDpeth_fast(cfg_peth, MUA, ipsi);  
    out_contra{iCell} = TSDpeth_fast(cfg_peth, MUA, contra);  
    
    peth_tvec(iCell,:) = out_ipsi{iCell}.tvec; X.peth_tvec{iCell} = peth_tvec(iCell,:);
    statindicestouse = peth_tvec(iCell,:) >= cfg_out.peth_Window(1) & peth_tvec(iCell,:) <= cfg_out.peth_Window(2);   % what window do we want to use for doing statistics on? The full peth window is just for display.
    statwindow = peth_tvec(iCell,statindicestouse);
    X.statwindow = statwindow;
    nums = find(statindicestouse);
    firstbin = min(nums); windowstart = peth_tvec(iCell,firstbin);
    lastbin = max(nums); windowend = peth_tvec(iCell,lastbin);
    binnum = length(out_ipsi{iCell}.tvec);
    binnumtouse = sum(statindicestouse); % sum up all the ones to get number of bins used for statistics
    stat_tvec = peth_tvec(iCell,statindicestouse);
    
    mean_FR_ipsi(iCell) = mean(out_ipsi{iCell}.data(statindicestouse)); 
    mean_FR_contra(iCell) = mean(out_contra{iCell}.data(statindicestouse)); 
    
    ipsi_Mean_peth = repmat(mean_FR_ipsi(iCell), 1, binnum);
    contra_Mean_peth = repmat(mean_FR_contra(iCell), 1, binnum);
    
    X.ipsi_Diff = abs(out_ipsi{iCell}.data(statindicestouse) - ipsi_Mean_peth(statindicestouse));
    X.contra_Diff = abs(out_contra{iCell}.data(statindicestouse) - contra_Mean_peth(statindicestouse));
    
    X.max_ipsi = max(X.ipsi_Diff);
    X.max_contra = max(X.contra_Diff);
    
    %% Caclulate the Circularly Shifted PETH
    tvec = tsdH.tvec; % tvec from the pupil position tsd
    ipsi_shuff{iCell} =  NaN(cfg_out.numShuff, binnumtouse);
    contra_shuff{iCell} = NaN(cfg_out.numShuff, binnumtouse);
    tic
    for iShuff = 1: cfg_out.numShuff
        h = waitbar(iShuff/cfg_out.numShuff);
        %     disp(iShuff);
        r(iShuff) = randsample(tvec,1); % choose a random value from the session times tvec
        ipsi_Shift = n + r(iShuff); % shift nasal times by a random amount
        contra_Shift = t + r(iShuff); % shift temporal times by a random amoun
        
        % SUBRACT SESSION DURATION FROM VALUES GREATER THAN THE LAST SESSION TIMESTAMP
        ipsi_List = ipsi_Shift > endtimetouse; ipsi_Indices = find(ipsi_List);
        contra_List = contra_Shift > endtimetouse; contra_Indices = find(contra_List); %
        
        ipsi_Shift_new = ipsi_Shift; contra_Shift_new = contra_Shift;
        if ~isempty(ipsi_Indices); ipsi_Shift_new = ipsi_Shift_new - endtimetouse; end  % subtract the session length from values > endtime
        if ~isempty(contra_Indices); contra_Shift_new = contra_Shift_new - endtimetouse; end  % subtract the session length from values > endtime
        
        % Calculate the shuffled PETH
        shuff_ipsi = TSDpeth_fast(cfg_peth, MUA, ipsi_Shift_new);
        ipsi_shuff{iCell}(iShuff,:) = shuff_ipsi.data(statindicestouse);
        shuff_ipsi_ave = nanmean(shuff_ipsi.data(statindicestouse));
        shuff_ipsi_uniform = repmat(shuff_ipsi_ave, 1, binnum);
        
        shuff_contra = TSDpeth_fast(cfg_peth, MUA, contra_Shift_new);
        contra_shuff{iCell}(iShuff,:) = shuff_contra.data(statindicestouse);
        shuff_contra_ave = nanmean(shuff_contra.data(statindicestouse));
        shuff_contra_uniform = repmat(shuff_contra_ave, 1, binnum);
        
        X.ipsi_Diff_shuff(iShuff,:) = abs(shuff_ipsi.data(statindicestouse) - shuff_ipsi_uniform(statindicestouse));  % *** not sure if I use true peth mean here or the shuffled peth mean.
        X.contra_Diff_shuff(iShuff,:) = abs(shuff_contra.data(statindicestouse) - shuff_contra_uniform(statindicestouse));
        
        X.max_ipsi_shuff(iShuff) = max(X.ipsi_Diff_shuff(iShuff,:));
        X.max_contra_shuff(iShuff) = max(X.contra_Diff_shuff(iShuff,:));
    end
    toc; disp('^^ time to calculate shuffled PETHs')
    
    set(groot,'ShowHiddenHandles','on')  % from matlab central
    delete(get(groot,'Children'))
    %% "Statistical test" - see where the sample value falls against the bootstrap distribution
    g_ipsi = X.max_ipsi_shuff < X.max_ipsi;  % logical of how many elements in the distribution are less than the max distance
    h_ipsi = nnz(g_ipsi)/length(X.max_ipsi_shuff);  % what is the fraction of distribution elements that are less than the max distance. e.g. 0.9 means the max distance is in the 90th percentile
    percent_ipsi(iCell) = h_ipsi;
    
    g_contra = X.max_contra_shuff < X.max_contra;  % logical of how many elements in the distribution are less than the max distance
    h_contra = nnz(g_contra)/length(X.max_contra_shuff);  % what is the fraction of distribution elements that are less than the max distance. e.g. 0.9 means the max distance is in the 90th percentile
    percent_contra(iCell) = h_contra;
    
    X.percent_ipsi = percent_ipsi(iCell);
    X.percent_contra = percent_contra(iCell);
    
    %% Determine whether the change in FR is an INCREASE or a DECREASE within the PETH window. 
    [max_ipsi, max_index_ipsi] = max(X.ipsi_Diff);
    bintouse = stat_tvec(max_index_ipsi);
    [~, index_ipsi] = min(abs(peth_tvec(iCell,:)-bintouse)); % find the bin in peth_tvec{iCell} which corresponds to the bin in stat_tvec that we know has the max difference
    ipsi_Peak = out_ipsi{iCell}.data(index_ipsi); 
    
    if ipsi_Peak > mean_FR_ipsi(iCell) 
        ipsi_Dir(iCell) = 1;  % Dir = direction of firing rate change. 1 is increased FR. -1 is decreased firing rate. 0 is equal. 
    elseif ipsi_Peak < mean_FR_ipsi(iCell)
        ipsi_Dir(iCell) = -1;
    elseif nPeak == mean_FR_ipsi(iCell)
        ipsi_Dir(iCell) = 0;
    else
        error('problem with determining firing rate change')
    end
    
    [max_contra, max_index_contra] = max(X.contra_Diff);
    bintouse = stat_tvec(max_index_contra);
    [~, index_contra] = min(abs(peth_tvec(iCell,:)-bintouse)); % find the bin in peth_tvec{iCell} which corresponds to the bin in stat_tvec that we know has the max difference
    contra_Peak = out_contra{iCell}.data(index_contra); 

    if contra_Peak > mean_FR_contra(iCell) 
        contra_Dir(iCell) = 1;  % Dir = direction of firing rate change. 1 is increased FR. -1 is decreased firing rate. 0 is equal. 
    elseif contra_Peak < mean_FR_contra(iCell)
        contra_Dir(iCell) = -1;
    elseif tPeak == mean_FR_contra(iCell)
        contra_Dir(iCell) = 0;
    else
        error('problem with determining firing rate change')
    end
    
    %% Plot It
    if cfg_out.doPlot == 1
        clf;
        subplot(321); % IPSI peth
        plot(peth_tvec(iCell,:), out_ipsi{iCell}.data);
        set(gca, 'FontSize', cfg_out.FontSize)
        ylabel('FR (Hz)')
        xlabel('time peri Saccade (s)'); c = axis;
        line([c(1) c(2)], [aveFR aveFR], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-'); c1 = axis;
        axis([cfg_out.window(1) cfg_out.window(2) 0 c1(4)])
        [~, cellID, ~] = fileparts(tfile{iCell});
        title(cellID)
        text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string', strcat('sess. ave FR = ', num2str(round(aveFR,1))))
        
        subplot(322) % temporal peth
        plot(peth_tvec(iCell,:), out_contra{iCell}.data, 'Color', 'r');
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
        legend('CONTRA', 'Location', 'Northwest')
        subplot(321);
        axis([cfg_out.window(1) cfg_out.window(2) 0 z])
        line([windowstart windowstart], [0 z], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', [.5 .5 .5])   % start of statistics window
        line([windowend windowend], [0 z], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', [.5 .5 .5])   % end of statistics window
        legend('IPSI', 'Location', 'Northwest')
        
        subplot(323); % nasal peth
        plot(peth_tvec(iCell,:), out_ipsi{iCell}.data);
        set(gca, 'FontSize', cfg_out.FontSize)
        ylabel('FR (Hz)')
        c3 = axis;
        axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) 0 c3(4)])
        line([c(1) c(2)], [mean_FR_ipsi(iCell) mean_FR_ipsi(iCell)], 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')   % aveFR line
        [a b] = max(X.ipsi_Diff);
        bintouse = stat_tvec(b);
        [val, index] = min(abs(peth_tvec(iCell,:)-bintouse)); % find the bin in peth_tvec{iCell} which corresponds to the bin in stat_tvec that we know has the max difference
        line([stat_tvec(b) stat_tvec(b)], [mean_FR_ipsi(iCell) out_ipsi{iCell}.data(index)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')   % plot a line showing the max diff btwn peth and aveFR
        
        
        subplot(324) % temporal peth
        plot(peth_tvec(iCell,:), out_contra{iCell}.data, 'Color', 'r');
        set(gca, 'FontSize', cfg_out.FontSize)
        c4 = axis;
        axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) 0 c4(4)])
        line([c(1) c(2)], [mean_FR_contra(iCell) mean_FR_contra(iCell)], 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')   % aveFR line
        [a b] = max(X.contra_Diff);
        bintouse = stat_tvec(b);
        [val, index] = min(abs(peth_tvec(iCell,:)-bintouse)); % find the bin in peth_tvec{iCell} which corresponds to the bin in stat_tvec that we know has the max difference
        line([stat_tvec(b) stat_tvec(b)], [mean_FR_contra(iCell) out_contra{iCell}.data(index)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')   % plot a line showing the max diff btwn peth and aveFR
        axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) 0 z])
        subplot(323);
        axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) 0 z])
        
        subplot(325) %
        hist(X.max_ipsi_shuff, cfg_out.histbinnum); c = axis;
        line([X.max_ipsi  X.max_ipsi], [c(3) c(4)], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '-')
        c5 = axis;
        xlabel('abs() max diff. of shuffles')
        ylabel('count')
        set(gca, 'FontSize', 15)
        title(strcat('percentile = ', num2str(percent_ipsi*100)))
        text('Units', 'Normalized', 'Position', [0.5, 0.8], 'string', strcat('num shuff = ', num2str(cfg_out.numShuff)))
        
        subplot(326) %
        hist(X.max_contra_shuff, cfg_out.histbinnum); c = axis;
        line([X.max_contra  X.max_contra], [c(3) c(4)], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '-')
        c6 = axis;
        xlabel('abs() max diff. of shuffles')
        ylabel('count')
        set(gca, 'FontSize', 15)
        z = max([c5(4) c6(4)]);
        axis([c6(1) c6(2) 0 z])
        line([X.max_contra  X.max_contra], [0 z], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '-')
        title(strcat('percentile = ', num2str(percent_contra*100)))
        text('Units', 'Normalized', 'Position', [0.5, 0.8], 'string', strcat('num shuff = ', num2str(cfg_out.numShuff)))
        
        subplot(325)
        axis([c5(1) c5(2) 0 z])
        line([X.max_ipsi  X.max_ipsi], [0 z], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '-')
%         pause(1)
%         close
    end
    % warning('on', 'all')
    disp(strcat('percent nasal = ', num2str(X.percent_ipsi)))
    disp(strcat('percent temporal = ', num2str(X.percent_contra)))
    Z{iCell} = X;
    clear X
end

nUP = ipsi_Dir == 1;                % cells with an increase in FR during nasal saccades, in the stats window 
nDOWN = ipsi_Dir == -1;             % cells with an decrease in FR during nasal saccades, in the stats window 

tUP = contra_Dir == 1;                % cells with an increase in FR during temporal saccades, in the stats window 
tDOWN = contra_Dir == -1;             % cells with an increase in FR during temporal saccades, in the stats window 

w = [];
%% ONE OR THE OTHER significant, sign of change not considered
w.contra_95 = percent_contra >= 0.95; w.contra_95_sum = sum(w.contra_95); w.contra_95_per = sum(w.contra_95)/numCells;     % temporal sig. @ 95th percentile 
w.contra_99 = percent_contra >= 0.99; w.contra_99_sum = sum(w.contra_99); w.contra_99_per = sum(w.contra_99)/numCells;     % temporal sig. @ 99th percentile  

w.ipsi_95 = percent_ipsi >= 0.95; w.ipsi_95_sum = sum(w.ipsi_95); w.ipsi_95_per = sum(w.ipsi_95)/numCells;                    % nasal sig. @ 95th percentile 
w.ipsi_99 = percent_ipsi >= 0.99; w.ipsi_99_sum = sum(w.ipsi_99); w.ipsi_99_per = sum(w.ipsi_99)/numCells;                    % nasal sig. @ 95th percentile 

w.contra_or_ipsi_95 = percent_contra >= 0.95 | percent_ipsi >= 0.95; w.sum_contra_or_ipsi_95 = sum(w.contra_or_ipsi_95);  w.per_contra_or_ipsi_95 = w.sum_contra_or_ipsi_95/numCells;
w.contra_or_ipsi_99 = percent_contra >= 0.99 | percent_ipsi >= 0.99; w.sum_contra_or_ipsi_99 = sum(w.contra_or_ipsi_99);  w.per_contra_or_ipsi_99 = w.sum_contra_or_ipsi_99/numCells;

%% BOTH  significant (one up and one down) 

% contra up and ipsi down _95
w.contra_up_and_ipsi_down_95 = percent_contra >= 0.95 & contra_Dir == 1 &  percent_contra >= 0.95 & ipsi_Dir == -1;   % Increasing nasal & decreasing temporal @ 95th percentile
w.sum_contra_up_and_ipsi_down_95 = sum(w.contra_up_and_ipsi_down_95); w.per_contra_up_and_ipsi_down_95 = w.sum_contra_up_and_ipsi_down_95/numCells;
% contra down and ipsi up_95
w.contra_down_and_ipsi_up_95 = percent_contra >= 0.95 & contra_Dir == -1 &  percent_contra >= 0.95 & ipsi_Dir == 1;   % vice versa
w.sum_contra_down_and_ipsi_up_95 = sum(w.contra_down_and_ipsi_up_95); w.per_contra_down_and_ipsi_up_95 = w.sum_contra_down_and_ipsi_up_95/numCells;
% contra up and ipsi down _99
w.contra_up_and_ipsi_down_99 = percent_contra >= 0.99 & contra_Dir == 1 &  percent_contra >= 0.99 & ipsi_Dir == -1;   % Increasing nasal & decreasing temporal @ 99th percentile
w.sum_contra_up_and_ipsi_down_99 = sum(w.contra_up_and_ipsi_down_99); w.per_contra_up_and_ipsi_down_99 = w.sum_contra_up_and_ipsi_down_99/numCells;
% contra down and ipsi up _99
w.contra_down_and_ipsi_up_99 = percent_contra >= 0.99 & contra_Dir == -1 &  percent_contra >= 0.99 & ipsi_Dir == 1;   % vice versa 
w.sum_contra_down_and_ipsi_up_99 = sum(w.contra_down_and_ipsi_up_99); w.per_contra_down_and_ipsi_up_99 = w.sum_contra_down_and_ipsi_up_99/numCells; 

%% Significant Increasing Firing Rate Response  (either saccade type, or both) 
w.contra_increase_95 = percent_contra >= 0.95 & contra_Dir == 1; w.sum_contra_increase_95 = sum(w.contra_increase_95); w.contra_increase_95_percentage = w.sum_contra_increase_95/numCells;
w.ipsi_increase_95 = percent_ipsi >= 0.95 & ipsi_Dir == 1; w.sum_ipsi_increase_95 = sum(w.ipsi_increase_95); w.ipsi_increase_95_percentage = w.sum_ipsi_increase_95/numCells;

w.contra_increase_99 = percent_contra >= 0.99 & contra_Dir == 1; w.sum_contra_increase_99 = sum(w.contra_increase_99); w.contra_increase_99_percentage = w.sum_contra_increase_99/numCells;
w.ipsi_increase_99 = percent_ipsi >= 0.99 & ipsi_Dir == 1; w.sum_ipsi_increase_99 = sum(w.ipsi_increase_99); w.ipsi_increase_99_percentage = w.sum_ipsi_increase_99/numCells;

f95_contra_incr = find(w.contra_increase_95);  % temporal increase - find indices
f99_contra_incr = find(w.contra_increase_99); 

f95_ipsi_incr = find(w.ipsi_increase_95);  % nasal increase - find indices
f99_ipsi_incr = find(w.ipsi_increase_99); 

both_incr_95 = intersect(f95_contra_incr, f95_ipsi_incr); num_both_incr_95 = length(both_incr_95); per_both_95 = num_both_incr_95/numCells;
both_incr_99 = intersect(f99_contra_incr, f99_ipsi_incr); num_both_incr_99 = length(both_incr_99); per_both_99 = num_both_incr_99/numCells;

%% Significant Decreasing Firing Rate Response   (either saccade type, or both) 
w.contra_decrease_95 = percent_contra >= 0.95 & contra_Dir == -1; w.sum_contra_decrease_95 = sum(w.contra_decrease_95); w.contra_decrease_95_percentage = w.sum_contra_decrease_95/numCells;
w.ipsi_decrease_95 = percent_ipsi >= 0.95 & ipsi_Dir == -1; w.sum_ipsi_decrease_95 = sum(w.ipsi_decrease_95); w.ipsi_decrease_95_percentage = w.sum_ipsi_decrease_95/numCells;

w.contra_decrease_99 = percent_contra >= 0.99 & contra_Dir == -1; w.sum_contra_decrease_99 = sum(w.contra_decrease_99); w.contra_decrease_99_percentage = w.sum_contra_decrease_99/numCells;
w.ipsi_decrease_99 = percent_ipsi >= 0.99 & ipsi_Dir == -1; w.sum_ipsi_decrease_99 = sum(w.ipsi_decrease_99); w.ipsi_decrease_99_percentage = w.sum_ipsi_decrease_99/numCells;

w.f95_contra_decr = find(w.contra_decrease_95);  % temporal increase - find indices
w.f99_contra_decr = find(w.contra_decrease_99); 

w.f95_ipsi_decr = find(w.ipsi_decrease_95);  % nasal increase - find indices
w.f99_ipsi_decr = find(w.ipsi_decrease_99); 

w.both_decr_95 = intersect(w.f95_contra_decr, w.f95_ipsi_decr); w.num_both_decr_95 = length(w.both_decr_95); w.per_both_decr_95 = w.num_both_decr_95/numCells;
w.both_decr_99 = intersect(w.f99_contra_decr, w.f99_ipsi_decr); w.num_both_decr_99 = length(w.both_decr_99); w.per_both_decr_99 = w.num_both_decr_99/numCells;

%% NEITHER   (not significant for either saccade type)
w.numSig_neither_95 = percent_contra < 0.95 & percent_ipsi < 0.95; w.sumSig_neither_95 = sum(w.numSig_neither_95); w.per_Sig_neither_95 = w.sumSig_neither_95/numCells;
w.numSig_neither_99 = percent_contra < 0.99 & percent_contra < 0.99; w.sumSig_neither_99 = sum(w.numSig_neither_99); w.per_Sig_neither_99 = w.sumSig_neither_99/numCells;

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