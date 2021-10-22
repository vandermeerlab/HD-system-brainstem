function [zVal_T, zVal_N, cellID, preSacWindwow, periSacWindow] = saccadePETHzscore(varargin)
% 4/2021. JJS.
% get Zscore value for nasal and temporal saccades for a given cell. Uses some -pre and -post saccade window for consideration.

periSacWindow = [-.05 0];
periDur = periSacWindow(2) - periSacWindow(1);
preSacWindwow = [-.2 -.1];
preDur = preSacWindwow(2) - preSacWindwow(1);
process_varargin(varargin)

cellCounterToUse = 0;
cellCounterIndex = 0;
fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN;
    disp(SSN)
    [S] = LoadSpikesJeff;
    if exist(strcat(SSN, '-VT1_proc.mat'))
        [temporalSaccades, nasalSaccades, ~, ~, ~, ~, ~, ~, ~] = processPupilData2([], 'doPlotEverything', 0, 'doPlotThresholds', 0); % Get saccade times
        %% Pre-Allocate
        numT = length(temporalSaccades);
        EventResponseT = nan(length(S.t),numT);
        meanEventResponseT = nan(1,length(S.t));
        
        numN = length(nasalSaccades);
        EventResponseN = nan(length(S.t),numN);
        meanEventResponseN = nan(1,length(S.t));
        for iCell = 1:length(S.t)
            cellCounterToUse = cellCounterToUse + 1;
            cellCounterIndex = cellCounterIndex +1;
            cellID(cellCounterToUse) = cellCounterIndex;
            %% Count the # of spikes during the event window for each TEMPROAL SACCADE
            for iT = 1:numT
                SrB = restrict(S, temporalSaccades(iT) + preSacWindwow(1), temporalSaccades(iT) + preSacWindwow(2));     % Spike train (S) restricted (r) to baseline (B) pre-saccade interval
                BaselineT(iCell,iT) = length(SrB.t{iCell}) / preDur;
                
                Sr = restrict(S, temporalSaccades(iT) + periSacWindow(1), temporalSaccades(iT) + periSacWindow(2));
                EventResponseT(iCell,iT) = length(Sr.t{iCell}) / periDur;
            end
            meanBaselineT(iCell) = mean(BaselineT(iCell,:));
            meanEventResponseT(iCell) = mean(EventResponseT(iCell,:));
            %% Count the # of spikes during the event window for each NASAL SACCADE
            for iT = 1:numN
                SrB = restrict(S, nasalSaccades(iT) + preSacWindwow(1), nasalSaccades(iT) + preSacWindwow(2));     % Spike train (S) restricted (r) to baseline (B) pre-saccade interval
                BaselineN(iCell,iT) = length(SrB.t{iCell}) / preDur;
                
                Sr = restrict(S, nasalSaccades(iT) + periSacWindow(1), nasalSaccades(iT) + periSacWindow(2));
                EventResponseN(iCell,iT) = length(Sr.t{iCell}) / periDur;
            end
            meanBaselineN(iCell) = mean(BaselineN(iCell,:));
            meanEventResponseN(iCell) = mean(EventResponseN(iCell,:));
            
            %% Calculate the zscores
            zVal_T(cellCounterToUse) = (meanEventResponseT(iCell) - mean(BaselineT(iCell,:)))/std(BaselineT(1,:)); %#ok<NODEF,*AGROW>
            zVal_N(cellCounterToUse) = (meanEventResponseN(iCell) - mean(BaselineN(iCell,:)))/std(BaselineN(1,:));
        end
    else
        disp('no eye tracking file detected. Skipping Session.')
        for iCell = 1:length(S.t)
            cellCounterIndex = cellCounterIndex +1;
        end
    end
    popdir;
end





