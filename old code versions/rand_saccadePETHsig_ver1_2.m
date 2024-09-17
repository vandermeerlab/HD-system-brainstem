function [out_n, out_t, mn_shuff, mt_shuff, peth_tvec, cfg_out] = rand_saccadePETHsig_ver1_2(tfile, cfg_in)
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
% ver1_2:   this version abandons SpikePETHvdm.m and instead uses tsdPETH with a previous call to MakeQfromS.


cfg_def.numShuff = 1000;
cfg_def.doPlot = 0;
cfg_out = ProcessConfig(cfg_def, cfg_in);

pushdir(fileparts(tfile));
temp = pwd; disp(temp)
%% Load Spikes
cfgS.uint = '64';
cfgS.fc = {tfile};
myCell = LoadSpikes(cfgS);
[csc_tsd, orientation, samplingrate, dt] = GetOrientationValues([]); %#ok<ASGLU>

%% Load Saccade times
load(FindFile('*saccades-edited.mat'), 'nasalSaccades', 'temporalSaccades', 'tsdH')
keep = find(~isnan(temporalSaccades)); t = temporalSaccades(keep);  %#ok<*FNDSB>
keep = find(~isnan(nasalSaccades)); n = nasalSaccades(keep);

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
cfg_peth.window = [-.2 .2];
cfg_peth.mode = 'interp';
out_n = TSDpeth_fast(cfg_peth, MUA_n, n);
out_t = TSDpeth_fast(cfg_peth, MUA_t, t);

%% Plot it
% if cfg_out.doPlot == 1
clf; plot(out_n); hold on; plot(out_t); pause(1); close 
% end

%% Caclulate the Random PETHs
tvec = tsdH.tvec; % tvec from the pupil position tsd
mt_shuff =  NaN(cfg_out.numShuff, length(out_n.tvec));
mn_shuff = NaN(cfg_out.numShuff, length(out_t.tvec));
tic
for iShuff = 1: cfg_out.numShuff
    waitbar(iShuff/cfg_out.numShuff)
%     disp(iShuff)
    % Generate the random times with length = num saccades
    tvecshort = tvec(length(out_t.tvec)+2: end - length(out_t.tvec) - 2); % shortening the tvec so that the PETH window does not run up against start or end of session (and have too few bins)
    for iN = 1:length(n)  % NASAL
        r_n(iN) = randsample(tvecshort,1);
    end
    for iT = 1:length(t) % TEMPORAL
        r_t(iT) = randsample(tvecshort,1);
    end
    %% Calculate the shuffled PETH
    temp_n = TSDpeth_fast(cfg_peth, MUA_n, r_n);
    temp_t = TSDpeth_fast(cfg_peth, MUA_t, r_t);
    
    mn_shuff(iShuff,:) = temp_n.data;
    mt_shuff(iShuff,:) = temp_t.data;
    
    peth_tvec = temp_n.tvec; % this should be the same as temp_t.tvec, and the same on every iteration of the for loop
end
toc; disp('^^ time to calculate shuffled PETHs')

if cfg_out.doPlot == 1
    clf; 
    plot(peth_tvec, nanmean(mn_shuff)); hold on; 
    plot(peth_tvec, nanmean(mt_shuff)); 
end

