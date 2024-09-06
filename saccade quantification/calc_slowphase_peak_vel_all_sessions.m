function [EVmov_filt_touse, maxEVfiltered, meanEVfiltered] = calc_slowphase_peak_vel_all_sessions(fd, varargin)

doSave = 1;
% doPlot = 1;
filtorder = 10;
if isempty(fd)
    fd = FindFiles('*keys.m');
end
disp(length(fd))
startSess = 1;
endSess = length(fd);
process_varargin(varargin);

maxEV = [];
for iSess = startSess: endSess
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN;
    filename = strcat(SSN, '-saccades-edited.mat');
        if exist(strcat(SSN, '-saccades-edited.mat')) == 2  % if the fig exists, then this session has already been done. *Hacky. Change to see if variable k exists.
            sd = LoadSessionData([]);
            [data_out, data_outR, nSpikesRemoved] = getSlowPhaseData([], sd, []);
            %% Remove Outlier Values
            EV = data_outR.horiz_eye_vel; 
            y = medfilt1(EV.data,filtorder,'omitnan');
            EVfilt = tsd(EV.tvec, y);
            %% Remove stationary periods
            load(strcat(SSN, '-AHV_StationaryTimes.mat'))
            EVmov_filt = antirestrict(EVfilt, STtstart, STtend);
%             zEVmov_filt_data = zscore(EVmov_filt.data); % regular z-score won't work if there are NaNs
            N = normalize(EVmov_filt.data); % z-score is the default method.  check to make sure mean = 0 and std = 1. M = mean(N, 'omitnan'); S = std(N, 'omitnan');           
            keepIndices = abs(N) < 2.5; % remove the tope 1% of the the (abs. values of the) distribution. 0.5% on each side. 
%             zEVmov_filtTSD = tsd(EVmov_filt.tvec(removeIndices), N(removeIndices)); % *** Add a plot as a check here. 
            EVmov_filt_touse{iSess} = tsd(EVmov_filt.tvec(keepIndices), EVmov_filt.data(keepIndices)); % plot(EVmov_filt_touse) 
            maxEVfiltered(iSess) = max(EVmov_filt_touse{iSess}.data);
            maxEVtouse = maxEVfiltered(iSess);
            disp(maxEVtouse); 
            meanEVfiltered(iSess) = nanmean(EVmov_filt_touse{iSess}.data); 
            meanEVtouse = meanEVfiltered(iSess); 
        else
            warning('saccades_edited.mat file not found')
        end
        if doSave
        save(filename, 'maxEV', 'maxEVtouse', 'meanEVtouse', "-append")
        end
end

