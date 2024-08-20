function [a,b, cfg_out] = saccadePETHsig_ver1_1(tfilelist, cfg_in)
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
a = []; b = [];
cfg_def.blank = [];
cfg_out = ProcessConfig(cfg_def, cfg_in);

for iC = 1:length(tfilelist)
    pushdir(fileparts(tfilelist{iC}));
    temp = pwd; disp(temp)
    %% Load Spikes
    cfgS.uint = '64';
    cfgS.fc = {tfilelist{iC}};
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
    
    starttime = csc_tsd.tvec(1);
    myCell.t{1} = myCell.t{1} - starttime;  % subtract the very first time stamp to convert from Unix time to 'start at zero' time.
    
    %% Load Saccade times
    load(FindFile('*saccades-edited.mat'), 'nasalSaccades', 'temporalSaccades')
    keep = find(~isnan(temporalSaccades)); t = temporalSaccades(keep);  %#ok<*FNDSB>
    keep = find(~isnan(nasalSaccades)); n = nasalSaccades(keep);
    
    cfg_peth = [];
    cfg_peth.dt = 0.01;
    cfg_peth.doPlot = 0;
    cfg_peth.window = [-.2 .2];
    [outputS_t, outputT_t, outputGau_t, outputIT_t, cfg_peth_t] = SpikePETHvdm(cfg_peth, myCell, t);  % temporal saccade peth
    [outputS_n, outputT_n, outputGau_n, outputIT_n, cfg_peth_n] = SpikePETHvdm(cfg_peth, myCell, n);  % nasal saccade peth
    
    mt = histcounts(outputS_t, outputIT_t);
    mn = histcounts(outputS_n, outputIT_n);
    
    plot((outputIT_t(1:end-1)), mt / cfg_peth.dt / length(t)); hold on    % *** hack. Using end-1 here because outputIT is 1 element longer. Using histc produced the right lenght, but the last entry was always zero. 
    plot(smoothdata(outputIT_n(1:end-1)), mn / cfg_peth.dt / length(n));
    
    
    
end
% m = histc(outputS, outputIT);
% bar(outputIT,m/cfg_peth.dt/length(t));
% 
% mN = histcounts(outputS, outputIT);
% FRxBinN(cellCounterToUse,:) = mN/cfg.dt/length(nasalSaccades);




% % generate a matrix of shuffled spike times
% h = waitbar(0, 'please wait');
% for iShuff = 1:numShuff
%     waitbar(iShuff/numShuff)
%     temp = ShuffleTS(cfg, myCell);
%     a.t{1} = temp.t{1};
%     a.cfg = myCell.cfg;
% %     shuffmatrix(iShuff,:) = ts_in.t{1}';
%     cfg_peth = [];
%     cfg_peth.dt = dt;
%     cfg_peth.doPlot = 0;
%     cfg_peth.window = [-.2 .2];
%     [outputS, outputT, outputGau, outputIT, cfg_peth] = SpikePETHvdm(cfg_peth, myCell, t);  % calculate spike peth on each shuffled spike train
%     M(iShuff,:) = histc(outputS, outputIT);  % arrange the results into a matrix with histc. These are spike counts.
%     Mfr(iShuff,:) = M(iShuff,:)/dt/length(t);   % these are firing rates
% %     bar(outputIT,m/cfg.dt/length(t));
%     if doPlot == 1
%         clf
% %         m = histc(outputS, outputIT);
%         bar(outputIT,m/cfg_peth.dt/length(t));
%     end
%     clear a
%     clear temp
% end