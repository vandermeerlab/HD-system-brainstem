function [out_n, out_t, peth_tvec, meanFR_n, meanFR_t, X, cfg_out] = ks_saccade_test_mean_ver1_1(tfile, cfg_in)
% saccadePETHsig_ver1_1.m  This function determines whether each neuorn is significantly modulated around the saccade time
%
%   Inputs
%           tfilelist: a cell array in which each element is the file location for a given neuron to use. For example,
%                        {'C:\Jeff\U01\datatouse\M282\M282-2022-02-04-1\M282-2022-02-04-1-TT01_2.t'}
%                       These will be all of the neurons to analyze. Most likely, the list of Confidence = 1 NPH neurons.
%           cfg_in: config variable.

%   Outputs
%           a:
%           b:
%
%   To do:
%           (1) write a function to go through all of the saccade files and remove the NaNs and re-save
% ver1_1    This function uses a ks test between the cumulative distributions of the actual peth and the mean FR of the peth.  

cfg_def.FontSize = 15;
cfg_def.numShuff = 1000; % how many shuffles to do 
cfg_def.doPlot = 1;
cfg_def.peth_Window = [-.5 .5]; % peth DISPLAY window
cfg_def.window = [-.21 .21]; % peth SIGNIFICANCE window 
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
statindicestouse = peth_tvec >= cfg_out.window(1) & peth_tvec <= cfg_out.window(2);   % what window do we want to use for doing statistics on? The full peth window is just for display. 
statwindow = peth_tvec(statindicestouse);
X.statwindow = statwindow;
nums = find(statindicestouse); 
firstbin = min(nums); windowstart = peth_tvec(firstbin);
lastbin = max(nums); windowend = peth_tvec(lastbin);
binnum = length(out_n.tvec);

meanFR_n = mean(out_n.data(statindicestouse));
meanFR_t = mean(out_t.data(statindicestouse));

nMean_peth = repmat(meanFR_n, 1, binnum);
tMean_peth = repmat(meanFR_t, 1, binnum);

%% Statistical test
[X.h_n, X.p_n, X.ks2stat_n] = kstest2(cumsum(out_n.data(statindicestouse)), cumsum(nMean_peth(statindicestouse)));
[X.h_t, X.p_t, X.ks2stat_t] = kstest2(cumsum(out_t.data(statindicestouse)), cumsum(tMean_peth(statindicestouse)));

if cfg_out.doPlot == 1
    clf;
%     f = figure;
    subplot(221); % nasal peth
%     f.Position = [600 400 945 700];  % does this not work with subplot?
    plot(peth_tvec, out_n.data); 
    text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string', strcat('ave FR = ', num2str(meanFR_n)))
    set(gca, 'FontSize', cfg_out.FontSize)
    ylabel('FR (Hz)')
    xlabel('time peri Saccade (s)')
    c = axis;
    line([windowstart windowend], [meanFR_n meanFR_n], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')   % aveFR line
    line([windowstart windowstart], [0 c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', [.5 .5 .5])   % start of statistics window
    line([windowend windowend], [0 c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', [.5 .5 .5])   % end of statistics window 
    axis([c(1) c(2) c(3) c(4)])
    [~, b, ~] = fileparts(tfile); 
    title(b)
    legend('nasal', 'Location', 'Northwest')
   
    subplot(222) % temporal peth
    plot(peth_tvec, out_t.data, 'Color', 'r'); 
    text('Units', 'Normalized', 'Position', [0.1, 0.1], 'string', strcat('ave FR = ', num2str(meanFR_t)))
    set(gca, 'FontSize', cfg_out.FontSize)
%     ylabel('FR (Hz)')
    xlabel('time peri Saccade (s)')
    c = axis;
    line([windowstart windowend], [meanFR_t meanFR_t], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-')
    line([windowstart windowstart], [0 c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', [.5 .5 .5])   % start of statistics window
    line([windowend windowend], [0 c(4)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--', 'Color', [.5 .5 .5])   % end of statistics window 
    title(b)
    legend('temporal', 'Location', 'Northwest')
    
    subplot(223) % nasal cumsum
    plot(peth_tvec(statindicestouse), cumsum(out_n.data(statindicestouse))); hold on;
    plot(peth_tvec(statindicestouse), cumsum(nMean_peth(statindicestouse)), 'k')
    legend('real', 'shuffle', 'Location', 'Northwest')
    set(gca, 'FontSize', cfg_out.FontSize)
    ylabel('cumsum FR (Hz)')
    c = axis; axis([windowstart windowend c(3) c(4)])
    xlabel(strcat('p = ', num2str(X.p_n)))
    
    subplot(224) % temporal cumsum
    plot(peth_tvec(statindicestouse), cumsum(out_t.data(statindicestouse)), 'r'); hold on;
    plot(peth_tvec(statindicestouse), cumsum(tMean_peth(statindicestouse)), 'k')
    legend('real', 'shuffle', 'Location', 'Northwest')
    set(gca, 'FontSize', cfg_out.FontSize)
    xlabel(strcat('p = ', num2str(X.p_t)))
    c = axis; axis([windowstart windowend c(3) c(4)])

end

