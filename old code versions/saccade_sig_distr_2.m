function [out_n, out_t, peth_tvec, meanFR_n, meanFR_t, mn_shuff, mt_shuff, X, cfg_out] = saccade_sig_distr_2(tfile, cfg_in)
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

cfg_def.FontSize = 15;
cfg_def.numShuff = 1000; % how many shuffles to do
cfg_def.histbinnum = 50;
cfg_def.doPlot = 1;
cfg_def.peth_Window = [-.2 .1]; % the window for doing statistics on. These values should be smaller than cfg_def.window
cfg_def.window = [-1 1]; % the window for display
cfg_out = ProcessConfig(cfg_def, cfg_in);

pushdir(fileparts(tfile));
temp = pwd; disp(temp)
%% Load Spikes
cfgS.uint = '64';
cfgS.fc = tfile;
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

cfg_peth = [];
cfg_peth.dt = 0.01;
cfg_peth.doPlot = 0;
cfg_peth.window = cfg_out.window; % this is the window for display purposes. the "peth_Window" is a smaller part of that for statistics.
cfg_peth.mode = 'interp';
out_n = TSDpeth_fast(cfg_peth, MUA_n, n);
out_t = TSDpeth_fast(cfg_peth, MUA_t, t);

peth_tvec = out_n.tvec; X.peth_tvec = peth_tvec;
statindicestouse = peth_tvec >= cfg_out.peth_Window(1) & peth_tvec <= cfg_out.peth_Window(2);   % what window do we want to use for doing statistics on? The full peth window is just for display.
statwindow = peth_tvec(statindicestouse);
X.statwindow = statwindow;
nums = find(statindicestouse);
firstbin = min(nums); windowstart = peth_tvec(firstbin);
lastbin = max(nums); windowend = peth_tvec(lastbin);
binnum = length(out_n.tvec);
binnumtouse = sum(statindicestouse); % sum up all the ones to get number of bins used for statistics
stat_tvec = peth_tvec(statindicestouse);

meanFR_n = mean(out_n.data(statindicestouse));
meanFR_t = mean(out_t.data(statindicestouse));

nMean_peth = repmat(meanFR_n, 1, binnum);
tMean_peth = repmat(meanFR_t, 1, binnum);

X.nDiff = abs(out_n.data(statindicestouse) - nMean_peth(statindicestouse));
X.tDiff = abs(out_t.data(statindicestouse) - tMean_peth(statindicestouse));

X.max_n = max(X.nDiff);
X.max_t = max(X.tDiff);

%% Caclulate the Circularly Shifted PETH
tvec = tsdH.tvec; % tvec from the pupil position tsd
mt_shuff =  NaN(cfg_out.numShuff, binnumtouse);
mn_shuff = NaN(cfg_out.numShuff, binnumtouse);
tic
for iShuff = 1: cfg_out.numShuff
    waitbar(iShuff/cfg_out.numShuff);
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
    mn_shuff(iShuff,:) = shuff_n.data(statindicestouse);
    shuff_n_ave = nanmean(shuff_n.data(statindicestouse));
    shuff_n_uniform = repmat(shuff_n_ave, 1, binnum);
    
    shuff_t = TSDpeth_fast(cfg_peth, MUA_t, tShift_new);
    mt_shuff(iShuff,:) = shuff_t.data(statindicestouse);
    shuff_t_ave = nanmean(shuff_t.data(statindicestouse));
    shuff_t_uniform = repmat(shuff_t_ave, 1, binnum);
    
    X.nDiff_shuff(iShuff,:) = abs(shuff_n.data(statindicestouse) - shuff_n_uniform(statindicestouse));  % *** not sure if I use true peth mean here or the shuffled peth mean.
    X.tDiff_shuff(iShuff,:) = abs(shuff_t.data(statindicestouse) - shuff_t_uniform(statindicestouse));
    
    X.max_n_shuff(iShuff) = max(X.nDiff_shuff(iShuff,:));
    X.max_t_shuff(iShuff) = max(X.tDiff_shuff(iShuff,:));
end
toc; disp('^^ time to calculate shuffled PETHs')
% close h
%% "Statistical test" - see where the sample value falls against the bootstrap distribution
g_n = X.max_n_shuff < X.max_n;  % logical of how many elements in the distribution are less than the max distance
h_n = nnz(g_n)/length(X.max_n_shuff);  % what is the fraction of distribution elements that are less than the max distance. e.g. 0.9 means the max distance is in the 90th percentile
percent_n = h_n;

g_t = X.max_t_shuff < X.max_t;  % logical of how many elements in the distribution are less than the max distance
h_t = nnz(g_t)/length(X.max_t_shuff);  % what is the fraction of distribution elements that are less than the max distance. e.g. 0.9 means the max distance is in the 90th percentile
percent_t = h_t;


%% Plot It
if cfg_out.doPlot == 1
    clf;
    subplot(321); % nasal peth
    plot(peth_tvec, out_n.data);
    set(gca, 'FontSize', cfg_out.FontSize)
    ylabel('FR (Hz)')
    xlabel('time peri Saccade (s)'); c = axis;
    line([c(1) c(2)], [aveFR aveFR], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-'); c1 = axis;
    axis([cfg_out.window(1) cfg_out.window(2) 0 c1(4)])
    [~, b, ~] = fileparts(tfile);
    title(b)
    text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string', strcat('sess. ave FR = ', num2str(round(aveFR,1))))
    
    subplot(322) % temporal peth
    plot(peth_tvec, out_t.data, 'Color', 'r');
    set(gca, 'FontSize', cfg_out.FontSize)
    xlabel('time peri Saccade (s)'); c = axis;
    line([c(1) c(2)], [aveFR aveFR], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-'); c2 = axis;
    axis([cfg_out.window(1) cfg_out.window(2) 0 c2(4)])
    title(b)
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
    plot(peth_tvec, out_n.data);
    set(gca, 'FontSize', cfg_out.FontSize)
    ylabel('FR (Hz)')
    c3 = axis;
    axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) 0 c3(4)])
    line([c(1) c(2)], [meanFR_n meanFR_n], 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')   % aveFR line
    [a b] = max(X.nDiff);
    bintouse = stat_tvec(b); 
    [val, index] = min(abs(peth_tvec-bintouse)); % find the bin in peth_tvec which corresponds to the bin in stat_tvec that we know has the max difference
    line([stat_tvec(b) stat_tvec(b)], [meanFR_n out_n.data(index)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')   % plot a line showing the max diff btwn peth and aveFR
    
    
    subplot(324) % temporal peth
    plot(peth_tvec, out_t.data, 'Color', 'r');
    set(gca, 'FontSize', cfg_out.FontSize)
    c4 = axis;
    axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) 0 c4(4)])
    line([c(1) c(2)], [meanFR_t meanFR_t], 'Color', [.5 .5 .5], 'LineWidth', 1, 'LineStyle', '--')   % aveFR line
    [a b] = max(X.tDiff);
    bintouse = stat_tvec(b); 
    [val, index] = min(abs(peth_tvec-bintouse)); % find the bin in peth_tvec which corresponds to the bin in stat_tvec that we know has the max difference
    line([stat_tvec(b) stat_tvec(b)], [meanFR_t out_t.data(index)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')   % plot a line showing the max diff btwn peth and aveFR
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
    
end

