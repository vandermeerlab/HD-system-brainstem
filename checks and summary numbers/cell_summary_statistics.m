function [confirmed_field, best_guess_field, X] = cell_summary_statistics(fd)
% JJS. 2024-08-20.
% This function tabulates the number of NPH and Gi neurons from all recording sessions in the input directory, as well as location confidence scores.

% Inputs:       fd - file directory. Usually this is a list of all recording sessions with eyetracking. Ifempty, then searches in current dir.

% Outputs:
%               confirmed_field  - 1 x nSess. Is ExpKeys.LesionStructureConfirmed present?  NaN = no result. 0 = not present. 1 = is present.
%               best_guess_field - 1 x nSess. Is ExpKeys.RecordingStructureBestGuess present?  NaN = no result. 0 = not present. 1 = is present.

%               X - structure with info for each brain structure
%
%                       .numNPH              -   # of NPH cells for each session (ranges from 0 - 4ish)
%                       .lesion_conf_NPH     -   confidence scores for each session, according to histolgoy evidence. Same as ExpKeys.LesionConfidence. 1 = best. 2 = moderate. 3 = some evidence, but uncertain.

%                       .numGuessNPH         -   # of NPH cells for each session w/ non-histology evidence. (ranges from 0 - 4ish)
%                       .guess_NPH           -   confidence score for each session, according to non-histology evidence. Same as ExpKeys.RecordingStructureBestGuess. 1 = best. 2 = moderate. 3 = some evidence, but uncertain.

if isempty(fd)
    fd = FindFiles('*keys.m');
end
confirmed_field = NaN(1,length(fd));  % pre-allocate
best_guess_field = NaN(1,length(fd));

for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; disp(SSN)
    if ~exist(strcat(SSN,'-VT1.mp4')); error('VIDEO FILE NOT FOUND'); end
    EvalKeys;
    
    confirmed_structureID = 0; 
    guess_structureID = 0; 
    %% ExpKeys.LesionStructureConfirmed
    % Look for matching strings                                             ***note: ANY DIFFERENCE IN SPELLING WILL RESULT IN A MISMATCH***
    if isfield(ExpKeys, 'LesionStructureConfirmed')
        confirmed_field(iSess) = 1;
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'NPH');               % nucleus prepositus
        lc_Gi = strcmp(ExpKeys.LesionStructureConfirmed, 'Gi');                 % gigantocellular nucleus               a.k.a. PGNRd in rats.
        lc_SGN = strcmp(ExpKeys.LesionStructureConfirmed, 'SGN');               % supragenual nucleus
        lc_MVN = strcmp(ExpKeys.LesionStructureConfirmed, 'MVN');               % medial vestibular nucleus
        lc_DTN = strcmp(ExpKeys.LesionStructureConfirmed, 'DTN');               % dorsal tegmental nucleus
        lc_Abducens = strcmp(ExpKeys.LesionStructureConfirmed, 'Abducens');     % abducens nucleus
        lc_6N = strcmp(ExpKeys.LesionStructureConfirmed, '6N');                 % abducens nucleus                      ***This is a repeat. Fix this.***
        lc_mlf = strcmp(ExpKeys.LesionStructureConfirmed, 'mlf');               % medial longitudinal fascicular        ***Change this to indicate something more precise.***
        lc_PDTg = strcmp(ExpKeys.LesionStructureConfirmed, 'PDTg');             % posterior dorsal tegmental nucleus    [***Maybe I should merge this with DTN?***]
        lc_LDTg = strcmp(ExpKeys.LesionStructureConfirmed, 'LDTg');             % lateral dorsal tegmental nucleus
        lc_CG = strcmp(ExpKeys.LesionStructureConfirmed, 'CG');                 % central grey
        lc_mRt = strcmp(ExpKeys.LesionStructureConfirmed, 'mRt');               % mesencephalic reticular formation
        lc_PNC = strcmp(ExpKeys.LesionStructureConfirmed, 'PNC');               % pontine reticular nucleus.            Gi turns into PNC in the anterior aspect
        
        % Sum it up
        X.numNPH_confirmed(iSess) = sum(lc_NPH);
        X.numGi_confirmed(iSess) = sum(lc_Gi);
        X.numSGN_confirmed(iSess) = sum(lc_SGN);
        X.numMVN_confirmed(iSess) = sum(lc_MVN);
        X.numDTN_confirmed(iSess) = sum(lc_DTN);
        X.numAbducens_confirmed(iSess) = sum(lc_Abducens);
        X.num6N_confirmed(iSess) = sum(lc_6N);
        X.nummlf_confirmed(iSess) = sum(lc_mlf);
        X.numPDTg_confirmed(iSess) = sum(lc_PDTg);
        X.numLDTg_confirmed(iSess) = sum(lc_LDTg);
        X.numCG_confirmed(iSess) = sum(lc_CG);
        X.nummRt_confirmed(iSess) = sum(lc_mRt);
        X.numPNC_confirmed(iSess) = sum(lc_PNC);
        X.numPNC_confirmed(iSess) = sum(lc_PNC);
    else
        warning('ExpKeys.LesionStructureConfirmed field is missing')
        confirmed_field(iSess) = 0;
        X.numNPH_confirmed(iSess) = NaN;
        X.numGi_confirmed(iSess) = NaN;
        X.numSGN_confirmed(iSess) = NaN;
        X.numMVN_confirmed(iSess) = NaN;
        X.numDTN_confirmed(iSess) = NaN;
        X.numAbducens_confirmed(iSess) = NaN;
        X.num6N_confirmed(iSess) = NaN;
        X.nummlf_confirmed(iSess) = NaN;
        X.numPDTg_confirmed(iSess) = NaN;
        X.numLDTg_confirmed(iSess) = NaN;
        X.numCG_confirmed(iSess) = NaN;
        X.nummRt_confirmed(iSess) = NaN;
        X.numPNC_confirmed(iSess) = NaN;
    end
    
    numConfirmed(iSess) = sum(X.numNPH_confirmed + X.numGi_confirmed + X.numSGN_confirmed + X.numMVN_confirmed + X.numDTN_confirmed + ...
        X.numAbducens_confirmed + X.num6N_confirmed + X.nummlf_confirmed + X.numPDTg_confirmed + X.numLDTg_confirmed + X.numCG_confirmed + ...
        X.nummRt_confirmed + X.numPNC_confirmed + X.numPNC_confirmed + X.numPNC_confirmed, 'omitnan');  % how many neurons present in the Confirmed field
    
    if numConfirmed(iSess) > 0; confirmed_structureID = 1; end  % 0 indicates no confirmed neurons. 1 indicates some confirmed neuron(s). 
        
    %% ExpKeys.RecordingStructureBestGuess
    % Look for matching strings                                                 ***note: ANY DIFFERENCE IN SPELLING WILL RESULT IN A MISMATCH***
    if isfield(ExpKeys, 'RecordingStructureBestGuess')
        best_guess_field(iSess) = 1;
        
        bg_NPH = strcmp(ExpKeys.RecordingStructureBestGuess, 'NPH');            % search for NPH string for best guess recordings
        bg_Gi = strcmp(ExpKeys.RecordingStructureBestGuess, 'Gi');              % search for Gi string for best guess recordings
        bg_SGN = strcmp(ExpKeys.RecordingStructureBestGuess, 'SGN');            % supragenual nucleus
        bg_MVN = strcmp(ExpKeys.RecordingStructureBestGuess, 'MVN');            % medial vestibular nucleus
        bg_DTN = strcmp(ExpKeys.RecordingStructureBestGuess, 'DTN');            % dorsal tegmental nucleus
        bg_Abducens = strcmp(ExpKeys.RecordingStructureBestGuess, 'Abducens');  % abducens nucleus
        bg_6N = strcmp(ExpKeys.RecordingStructureBestGuess, '6N');              % abducens nucleus                      ***This is a repeat. Fix this.
        bg_mlf = strcmp(ExpKeys.RecordingStructureBestGuess, 'mlf');            % medial longitudinal fascicular        ***Change this to indicate something more precise.
        bg_PDTg = strcmp(ExpKeys.RecordingStructureBestGuess, 'PDTg');          % posterior dorsal tegmental nucleus    [***Maybe I should merge this with DTN?]
        bg_LDTg = strcmp(ExpKeys.RecordingStructureBestGuess, 'LDTg');          % lateral dorsal tegmental nucleus
        bg_CG = strcmp(ExpKeys.RecordingStructureBestGuess, 'CG');              % central grey
        bg_mRt = strcmp(ExpKeys.RecordingStructureBestGuess, 'mRt');            % mesencephalic reticular formation
        bg_PNC = strcmp(ExpKeys.RecordingStructureBestGuess, 'PNC');            % pontine reticular nucleus
        bg_unknown = strcmp(ExpKeys.RecordingStructureBestGuess, 'unknown');    % no good information is present 

        % Sum it up
        X.numNPH_guess(iSess) = sum(bg_NPH);
        X.numGi_guess(iSess) = sum(bg_Gi);
        X.numSGN_guess(iSess) = sum(bg_SGN);
        X.numMVN_guess(iSess) = sum(bg_MVN);
        X.numDTN_guess(iSess) = sum(bg_DTN);
        X.numAbducens_guess(iSess) = sum(bg_Abducens);
        X.num6N_guess(iSess) = sum(bg_6N);
        X.nummlf_guess(iSess) = sum(bg_mlf);
        X.numPDTg_guess(iSess) = sum(bg_PDTg);
        X.numLDTg_guess(iSess) = sum(bg_LDTg);
        X.numCG_guess(iSess) = sum(bg_CG);
        X.nummRt_guess(iSess) = sum(bg_mRt);
        X.numPNC_guess(iSess) = sum(bg_PNC);
        X.numunknown_guess(iSess) = sum(bg_unknown); 
    else
        warning('ExpKeys.RecordingStructureBestGuess field is absent.')
        best_guess_field(iSess) = 0;
        X.numNPH_guess(iSess) = NaN;
        X.numGi_guess(iSess) = NaN;
        X.numSGN_guess(iSess) = NaN;
        X.numMVN_guess(iSess) = NaN;
        X.numDTN_guess(iSess) = NaN;
        X.numAbducens_guess(iSess) = NaN;
        X.num6N_guess(iSess) = NaN;
        X.nummlf_guess(iSess) = NaN;
        X.numPDTg_guess(iSess) = NaN;
        X.numLDTg_guess(iSess) = NaN;
        X.numCG_guess(iSess) = NaN;
        X.nummRt_guess(iSess) = NaN;
        X.numPNC_guess(iSess) = NaN;
        X.numunknown_guess(iSess) = NaN; 
    end
    
        numGuess(iSess) = sum(X.numNPH_guess + X.numGi_guess + X.numSGN_guess + X.numMVN_guess + X.numDTN_guess + ...
        X.numAbducens_guess + X.num6N_guess + X.nummlf_guess + X.numPDTg_guess + X.numLDTg_guess + X.numCG_guess + ...
        X.nummRt_guess + X.numPNC_guess + X.numPNC_guess + X.numPNC_guess, 'omitnan');  % how many neurons present in the BestGuess field
    
        if numGuess(iSess) > 0; guess_structureID = 1; end   % 0 indicates no confirmed neurons. 1 indicates some confirmed neuron(s). 
        
        if confirmed_structureID + guess_structureID == 1
            warning('problem with structureID keys fields');  % One of these fields should have a strucutre ID for one or more neurons. Both fields should NOT have entries in the same keys file. 
        end
    clear ExpKeys
    popdir;
end

%% Calculate the sum total over sessions -- CONFIRMED
X.sumNPH_confirmed = nansum(X.numNPH_confirmed);
X.sumGi_confirmed = nansum(X.numGi_confirmed);
X.sumSGN_confirmed = nansum(X.numSGN_confirmed);
X.sumMVN_confirmed = nansum(X.numMVN_confirmed);
X.sumDTN_confirmed = nansum(X.numDTN_confirmed);
X.sumAbducens_confirmed = nansum(X.numAbducens_confirmed);
X.sum6N_confirmed = nansum(X.num6N_confirmed);
X.summlf_confirmed = nansum(X.nummlf_confirmed);
X.sumPDTg_confirmed = nansum(X.numPDTg_confirmed);
X.sumLDTg_confirmed = nansum(X.numLDTg_confirmed);
X.sumCG_confirmed = nansum(X.numCG_confirmed);
X.summRt_confirmed = nansum(X.nummRt_confirmed);
X.sumPNC_confirmed = nansum(X.numPNC_confirmed);

%% Calculate the sum total over sessions -- BEST GUESS
X.sumNPH_guess = nansum(X.numNPH_guess);
X.sumGi_guess = nansum(X.numGi_guess);
X.sumSGN_guess = nansum(X.numSGN_guess);
X.sumMVN_guess = nansum(X.numMVN_guess);
X.sumDTN_guess = nansum(X.numDTN_guess);
X.sumAbducens_guess = nansum(X.numAbducens_guess);
X.sum6N_guess = nansum(X.num6N_guess);
X.summlf_guess = nansum(X.nummlf_guess);
X.sumPDTg_guess = nansum(X.numPDTg_guess);
X.sumLDTg_guess = nansum(X.numLDTg_guess);
X.sumCG_guess = nansum(X.numCG_guess);
X.summRt_guess = nansum(X.nummRt_guess);
X.sumPNC_guess = nansum(X.numPNC_guess);
X.sumunknown_guess = nansum( X.numunknown_guess);

%% Calculate the total over both categories -- CONFIRMED + BESTGUESS 
X.NPH = X.sumNPH_confirmed + X.sumNPH_guess;
X.Gi = X.sumGi_confirmed + X.sumGi_guess;
X.SGN = X.sumSGN_confirmed + X.sumSGN_guess;
X.MVN = X.sumMVN_confirmed + X.sumMVN_guess;
X.DTN = X.sumDTN_confirmed + X.sumDTN_guess;
X.Abducens = X.sumAbducens_confirmed + X.sumAbducens_guess + X.sum6N_confirmed + X.sum6N_guess;
X.mlf = X.summlf_confirmed + X.summlf_guess;
X.PDTg = X.sumPDTg_confirmed + X.sumPDTg_guess;
X.LDTg = X.sumLDTg_confirmed + X.sumLDTg_guess;
X.CG = X.sumCG_confirmed + X.sumCG_guess;
X.mRt = X.summRt_confirmed + X.summRt_guess;
X.PNC = X.sumPNC_confirmed + X.sumPNC_guess;
X.unknown = X.sumunknown_guess;

X.total_neurons = X.NPH + X.Gi + X.SGN + X.MVN + X.DTN + X.Abducens + X.mlf + X.PDTg + X.LDTg + X.CG + X.mRt + X.PNC + X.unknown; 

























%% Check for Confidence Score from lesion evidence

%     if isfield(ExpKeys, 'LesionStructureConfirmed')
%         if isfield(ExpKeys, 'LesionConfidence')        % see if 'LesionConfidence' exists in ExpKeys.  ***It should in every session.
%             if isempty(ExpKeys.LesionConfidence)  % This field shoule be a cell {}. See if it is empty.
%                 lesion_conf_NPH(iSess) = NaN;  % if empty, then ExpKeys.RecrodingStructureBestGuess should be NOT empty.
%             else
%                 lesion_conf_NPH(iSess) = ExpKeys.LesionConfidence;
%             end
%         else
%             lesion_conf_NPH(iSess) = NaN; disp('ExpKeys.LesionConfidence field is not present.')
%         end
%     else
%         disp('ExpKeys.LesionConfidence field is not present.')
%     end


%         if isfield(ExpKeys, 'RecordingStructureBestGuess')
%             if isfield(ExpKeys, 'NeuronConfidence')        % see if 'LesionConfidence' exists in ExpKeys.  ***It should in every session.
%                 if isempty(ExpKeys.NeuronConfidence)  % This field shoule be a cell {}. See if it is empty.
%                     guess_conf_NPH(iSess) = NaN;  % if empty, then ExpKeys.LesionStructureConfirmed should be NOT empty.
%                 else
%                     guess_conf_NPH(iSess) = ExpKeys.NeuronConfidence;
%                 end
%             else
%                 guess_conf_NPH(iSess) = NaN; disp('ExpKeys.NeuronConfidence field is not present.')
%             end
%         else
%             disp('ExpKeys.RecordingStructureBestGuess field is not present.')
%         end
%
%         %% Video is absent. Enter NaNs.
%     else
%         numNPH_confirmed(iSess) = NaN;
%         numGi_confirmed(iSess) = NaN;
%         lesion_conf_NPH(iSess) = NaN;
%         guess_conf_NPH(iSess) = NaN;
%         disp('No eyetracking video found for this session.')
%     end
%
%     if isnan(lesion_conf_NPH(iSess)) && isnan(guess_conf_NPH(iSess))
%         warning('NaNs: THIS SESSION HAS NEITHER A LESION CONFIDENCE SCORE NOR A BEST GUESS CONFIDENCE SCORE')
%     end
%     if mod(lesion_conf_NPH(iSess),1)==1 && mod(guess_conf_NPH(iSess),1) == 1   % why doesn't this one work?
%         warning('mod ~= 0: THIS SESSION HAS NEITHER A LESION CONFIDENCE SCORE NOR A BEST GUESS CONFIDENCE SCORE')
%     end















% Z.numNPH_confirmed = numNPH_confirmed; Z.sumNPH_confirmed = sum(Z.numNPH_confirmed);
% Z.numGi_confirmed = numGi_confirmed;   Z.sumGi_confirmed = sum(Z.numGi_confirmed);
%
% Z.numNPH_best_guess = numNPH_best_guess; Z.sumNPH_best_guess = sum(Z.numNPH_best_guess);
% Z.numGi_best_guess = numGi_best_guess;   Z.sumGi_best_guess = sum(Z.numGi_best_guess);


% %% Tabulate the Results for NPH
% Z.lesion_conf_NPH = lesion_conf_NPH;
%
% % Level 1 confidence
% Z.lesion_conf_NPH_1 = lesion_conf_NPH == 1; % Length 1 x n Sessions. Confidence value for each session.
% Z.numNPH_confirmed = numNPH_confirmed(Z.lesion_conf_NPH_1); % Length 1 x n Sessions. # of NPH neurons recorded in each session.
% Z.NPH_1 = sum(Z.numNPH_confirmed); % # of conf. 1 NPH neurons.
%
% % Level 2 confidence
% Z.lesion_conf_NPH_2 = lesion_conf_NPH == 2; % Length 1 x n Sessions. Confidence value for each session.
% Z.numNPH_confirmed_2 = numNPH_confirmed(Z.lesion_conf_NPH_2); % Length 1 x n Sessions. # of NPH neurons recorded in each session.
% Z.NPH_2 = sum(Z.numNPH_confirmed_2); % # of conf. 1 NPH neurons.
%
% % Level 3 confidence
% Z.lesion_conf_NPH_3 = lesion_conf_NPH == 3; % Length 1 x n Sessions. Confidence value for each session.
% Z.numNPH_confirmed_3 = numNPH_confirmed(Z.lesion_conf_NPH_3); % Length 1 x n Sessions. # of NPH neurons recorded in each session.
% Z.NPH_3 = sum(Z.numNPH_confirmed_3); % # of conf. 1 NPH neurons.
%
% %% Tabulate the Results for Gi
% Z.lesion_conf_Gi = lesion_conf_Gi;
%
% % Level 1 confidence
% Z.lesion_conf_Gi_1 = lesion_conf_Gi == 1; % Length 1 x n Sessions. Confidence value for each session.
% Z.numGi_confirmed = numGi_confirmed(Z.lesion_conf_Gi_1); % Length 1 x n Sessions. # of NPH neurons recorded in each session.
% Z.NPH_1 = sum(Z.numNPH_confirmed); % # of conf. 1 NPH neurons.
%
% % Level 2 confidence
% Z.lesion_conf_Gi_2 = lesion_conf_Gi == 2; % Length 1 x n Sessions. Confidence value for each session.
% Z.numGi_confirmed_2 = numGi_confirmed(Z.lesion_conf_Gi_2); % Length 1 x n Sessions. # of NPH neurons recorded in each session.
% Z.Gi_2 = sum(Z.numGi_confirmed_2); % # of conf. 1 NPH neurons.
%
% % Level 3 confidence
% Z.lesion_conf_Gi_3 = lesion_conf_Gi == 3; % Length 1 x n Sessions. Confidence value for each session.
% Z.numGi_confirmed_3 = numGi_confirmed(Z.lesion_conf_Gi_3); % Length 1 x n Sessions. # of NPH neurons recorded in each session.
% Z.Gi_3 = sum(Z.numGi_confirmed_3); % # of conf. 1 NPH neurons.



%% Print the results











