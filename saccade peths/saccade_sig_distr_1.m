function [out_n, out_t, peth_tvec, meanFR_n, meanFR_t, X, cfg_out] = saccade_sig_distr_1(tfile, cfg_in)
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

cfg_def.FontSize = 15;
cfg_def.numShuff = 1000; % how many shuffles to do
cfg_def.doPlot = 1;
cfg_def.peth_Window = [-.2 .05];
% cfg_def.window = [-.21 .21];
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
cfg_peth.window = cfg_out.peth_Window;
cfg_peth.mode = 'interp';
out_n = TSDpeth_fast(cfg_peth, MUA_n, n);
out_t = TSDpeth_fast(cfg_peth, MUA_t, t);

peth_tvec = out_n.tvec;
% statindicestouse = peth_tvec >= cfg_out.peth_Window(1) & peth_tvec <= cfg_out.peth_Window(2);   % what window do we want to use for doing statistics on? The full peth window is just for display.
% statwindow = peth_tvec(statindicestouse);
% X.statwindow = statwindow;
% nums = find(statindicestouse);
% firstbin = min(nums); windowstart = peth_tvec(firstbin);
% lastbin = max(nums); windowend = peth_tvec(lastbin);
binnum = length(out_n.tvec);

meanFR_n = mean(out_n.data);
meanFR_t = mean(out_t.data);

nMean_peth = repmat(meanFR_n, 1, binnum);
tMean_peth = repmat(meanFR_t, 1, binnum);

X.nDiff = out_n.data - nMean_peth;
X.tDiff = out_t.data - tMean_peth;

X.max_n = max(X.nDiff);
X.max_t = max(X.tDiff);

%% Caclulate the Circularly Shifted PETH
tvec = tsdH.tvec; % tvec from the pupil position tsd
mt_shuff =  NaN(cfg_out.numShuff, length(out_n.tvec));
mn_shuff = NaN(cfg_out.numShuff, length(out_t.tvec));
tic
for iShuff = 1: cfg_out.numShuff
    waitbar(iShuff/cfg_out.numShuff)
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
    mn_shuff(iShuff,:) = shuff_n.data;
    shuff_n_ave = nanmean(shuff_n.data);
    shuff_n_uniform = repmat(shuff_n_ave, 1, binnum);
    
    shuff_t = TSDpeth_fast(cfg_peth, MUA_t, tShift_new);
    mt_shuff(iShuff,:) = shuff_t.data;
    shuff_t_ave = nanmean(shuff_t.data);
    shuff_t_uniform = repmat(shuff_t_ave, 1, binnum);
    
    X.nDiff_shuff(iShuff,:) = shuff_n.data - shuff_n_uniform;  % *** not sure if I use true peth mean here or the shuffled peth mean.
    X.tDiff_shuff(iShuff,:) = shuff_t.data - shuff_t_uniform;
    
    X.max_n_shuff(iShuff) = max(abs(X.nDiff_shuff(iShuff,:)));
    X.max_t_shuff(iShuff) = max(abs(X.tDiff_shuff(iShuff,:)));
end
toc; disp('^^ time to calculate shuffled PETHs')

%% "Statistical test" - see where the sample value falls against the bootstrap distribution

%% Plot It
if cfg_out.doPlot == 1
    clf;
    subplot(221); % nasal peth
    plot(peth_tvec, out_n.data);
    text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string', strcat('sess. ave FR = ', num2str(aveFR)))
    set(gca, 'FontSize', cfg_out.FontSize)
    ylabel('FR (Hz)')
    xlabel('time peri Saccade (s)')
    c = axis;
    axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) c(3) c(4)])
    line([c(1) c(2)], [aveFR aveFR], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')   % aveFR line
    [~, b, ~] = fileparts(tfile);
    title(b)
    legend('nasal', 'Location', 'Northwest')
    
    subplot(222) % temporal peth
    plot(peth_tvec, out_t.data, 'Color', 'r');
    text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string', strcat('sess. ave FR = ', num2str(aveFR)))
    set(gca, 'FontSize', cfg_out.FontSize)
    xlabel('time peri Saccade (s)')
    c = axis;
    axis([cfg_out.peth_Window(1) cfg_out.peth_Window(2) c(3) c(4)])
    line([c(1) c(2)], [aveFR aveFR], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')   % aveFR line
    title(b)
    legend('temporal', 'Location', 'Northwest')
    
    subplot(223) %
    hist(X.max_n_shuff, 20)
    c = axis;
    line([X.max_n  X.max_n], [c(3) c(4)], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '-')   
    xlabel('max diff. of shuffles')
    ylabel('count')
    set(gca, 'FontSize', 15)
    
    subplot(224) %
    hist(X.max_t_shuff, 20)
    c = axis;
    line([X.max_t  X.max_t], [c(3) c(4)], 'Color', 'g', 'LineWidth', 1, 'LineStyle', '-')
    xlabel('max diff. of shuffles')
    ylabel('count')
    set(gca, 'FontSize', 15)
end

