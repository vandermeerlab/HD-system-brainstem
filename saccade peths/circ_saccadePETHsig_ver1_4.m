function [c, tNew, nNew, mt_shuff, mn_shuff, outputS_t_shuff, outputS_n_shuff, outputIT_t_shuff, outputIT_n_shuff, cfg_out] = circ_saccadePETHsig_ver1_4(tfilelist, cfg_in)
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
% ver1_3:   combined the multiple for loops over shuffles into one for loop.
% ver1_4:   here I did circular shift of the order of the saccade indices. This is incorrect and accomplishes nothing. 

cfg_def.numShuff = 1000;
cfg_def.doPlot = 1;
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
    load(FindFile('*saccades-edited.mat'), 'nasalSaccades', 'temporalSaccades', 'tsdH')
    keep = find(~isnan(temporalSaccades)); t = temporalSaccades(keep);  %#ok<*FNDSB>
    keep = find(~isnan(nasalSaccades)); n = nasalSaccades(keep);
    
    %% Calculate the true PETH
    cfg_peth = [];
    cfg_peth.dt = 0.01;
    cfg_peth.doPlot = 0;
    cfg_peth.window = [-.2 .2];
    [outputS_t, ~, ~, outputIT_t, cfg_peth_t] = SpikePETHvdm(cfg_peth, myCell, t);  % temporal saccade peth
    [outputS_n, ~, ~, outputIT_n, cfg_peth_n] = SpikePETHvdm(cfg_peth, myCell, n);  % nasal saccade peth
    
    mt = histcounts(outputS_t, outputIT_t);
    mn = histcounts(outputS_n, outputIT_n);
    
    %% Plot it
    if cfg_out.doPlot == 1
        clf
        plot(outputIT_t(1:end-1), mt / cfg_peth.dt / length(t)); hold on    % *** hack. Using end-1 here because outputIT is 1 element longer. Using histc produced the right lenght, but the last entry was always zero.
        plot(outputIT_n(1:end-1), mn / cfg_peth.dt / length(n));
        %         disp('press any key to continue'); pause
        close
    end
    mt_shuff =  NaN(cfg_out.numShuff, length(outputIT_t)-1);
    mn_shuff = NaN(cfg_out.numShuff, length(outputIT_t)-1);
    %% Find the tvec indices for the saccades
    tvec = tsdH.tvec; % tvec from the pupil position tsd
    for iN = 1:length(n)  % NASAL
        closest = interp1(tvec, tvec, n(iN), 'nearest');   % find closest tvec value to this specific saccade
        indextouse_N(iN) = find(tvec == closest);
    end
    for iT = 1:length(t) % TEMPORAL
        closest = interp1(tvec, tvec, t(iT), 'nearest');   % find closest tvec value to this specific saccade
        indextouse_T(iT) = find(tvec == closest);
    end
    %% Perform cicular shift on the SACCADE index
    tic
    c_n = NaN(cfg_out.numShuff, length(indextouse_N));
    c_t = NaN(cfg_out.numShuff, length(indextouse_T));
    circ_n_timestamps = NaN(cfg_out.numShuff, length(indextouse_N));
    circ_t_timestamps = NaN(cfg_out.numShuff, length(indextouse_T));
    for iShuff = 1: cfg_out.numShuff
        %         disp(iShuff)
        waitbar(iShuff/cfg_out.numShuff)
        r_n = randsample(1:indextouse_N,1);
        r_t = randsample(1:indextouse_N,1);
        c_n(iShuff,:) = circshift(indextouse_N, r_n); % nasal
        c_t(iShuff,:) = circshift(indextouse_T, r_t); % temporal
        %% Convert back to timestamps
        circ_n_timestamps(iShuff,:) =  tvec(c_n(iShuff,:)); % nasal
        circ_t_timestamps(iShuff,:) =  tvec(c_t(iShuff,:)); % temporal
        
        %% Calculate the shuffled PETH
        [outputS_n_shuff, ~, ~, outputIT_n_shuff, ~] = SpikePETHvdm(cfg_peth, myCell, circ_n_timestamps(iShuff,:));  % nasal saccade peth
        [outputS_t_shuff, ~, ~, outputIT_t_shuff, ~] = SpikePETHvdm(cfg_peth, myCell, circ_t_timestamps(iShuff,:));  % temporal saccade peth
        
        mn_shuff(iShuff,:) = histcounts(outputS_n_shuff, outputIT_n_shuff); % nasal
        mt_shuff(iShuff,:) = histcounts(outputS_t_shuff, outputIT_t_shuff); % temporal
    end
    toc; disp('^^ time to calculate shuffled PETHs')
    
    if cfg_out.doPlot == 1
        figure
        plot(outputIT_t_shuff(1:end-1), mt_shuff / cfg_peth.dt / length(t)); hold on    % *** hack. Using end-1 here because outputIT is 1 element longer. Using histc produced the right lenght, but the last entry was always zero.
        plot(outputIT_n_shuff(1:end-1), mn_shuff / cfg_peth.dt / length(n));
        close
    end
end
