function [zVal_T, zVal_N, sig, cellID, preSacWindwow, periSacWindow] = saccadePETHzscore(cellID, varargin)
% 4/2021. JJS.
% get Zscore value for nasal and temporal saccades for a given cell. Uses some -pre and -post saccade window for consideration.
% cellID = the order of each cell considered here (those with video data) with respect to the set of all cells (a larger group). 

periSacWindow = [-.05 .05];
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

% calc. sig. cells
sig.ind_N = abs(zVal_N_old)>2;
sig.ind_T = abs(zVal_T_old)>2; 
sig.num_N = sum(sig.ind_N);
sig.num_T = sum(sig.ind_T);
sig.either = sig.ind_N | sig.ind_T;
sig.num_either = sum(sig.ind_N | sig.ind_T); 

[~, saccadeSigIDs] = find(sig.either); 
saccadeSigIDall = cellID(saccadeSigIDs);
%% 
% cellID = neurons that come from sessions with video tracking data
cellIDlogical = zeros(1,length(numCells));
cellIDlogical(cellID) = 1;
qsi_all = cellID(qual_sig_index); 
QSI = zeros(1,length(numCells));
QSI(qsi_all) = 1;

overlap_Asym = QSI & Y.A; 
overlap_Sym = QSI & (Y.V | Y.C);
overlap = QSI & Y.S; 

Asymmetric_w_video = Y.A & cellIDlogical; 
Symmetric_w_video = (Y.V | Y.C) & cellIDlogical; 


Asymmetric_and_Saccade = (QSI & Y.S); 
percent_Asym_overlap = sum(Asymmetric_and_Saccade)/sum(Asymmetric_w_video); 





