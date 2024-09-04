function [Z] = cell_summary_statistics(fd)
% JJS. 2024-08-20.
% This function tabulates the number of NPH and Gi neurons from all recording sessions in the input directory, as well as location confidence scores.

% Inputs:       fd - file directory. Usually this is a list of all recording sessions with eyetracking. Ifempty, then searches in current dir.

% Outputs:
%
%           numNPH              -   # of NPH cells for each session (ranges from 0 - 4ish)
%           numGi               -   # of Gi cells for each session
%           lesion_conf_NPH     -   confidence scores for each session, according to histolgoy evidence. Same as ExpKeys.LesionConfidence. 1 = best. 2 = moderate. 3 = some evidence, but uncertain.
%           guess_NPH           -   confidence score for each session, according to non-histology evidence. Same as ExpKeys.RecordingStructureBestGuess. 1 = best. 2 = moderate. 3 = some evidence, but uncertain.
%           numGuessNPH         -   # of NPH cells for each session w/ non-histology evidence.
%           NPHguess_conf       -
%           Z.
%
%
%
%

if isempty(fd)
    fd = FindFiles('*keys.m');
end
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; disp(SSN)
    if exist(strcat(SSN,'-VT1.mp4')) 
%         disp('video file found')
        EvalKeys;
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'NPH');       % nucleus prepositus
        lc_Gi = strcmp(ExpKeys.LesionStructureConfirmed, 'Gi');         % gigantocellular nucleus               a.k.a. PGNRd in rats. 
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'SGN');       % gupragenual nucleus
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'MVN');       % medial vestibular nucleus
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'DTN');       % dorsal tegmental nucleus
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'Abducens');  % abducens nucleus
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, '6N');        % abducens nucleus                      ***This is a repeat. Fix this.  
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'mlf');       % medial longitudinal fascicular        ***Change this to indicate something more precise.
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'PDTg');      % posterior dorsal tegmental nucleus    [***Maybe I should merge this with DTN?]
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'LDTg');      % lateral dorsal tegmental nucleus 
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'CG');        % central grey
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'mRt');       % mesencephalic reticular formation
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'PNC');       % pontine reticular nucleus.            Gi turns into PNC in the anterior aspect 
   
        
        
        
        
        bg_NPH = strcmp(ExpKeys.RecordingStructureBestGuess, 'NPH');   % search for NPH string for best guess recordings
        bg_Gi = strcmp(ExpKeys.RecordingStructureBestGuess, 'Gi');     % search for Gi string for best guess recordings
        
        
        
        
        numNPH_confirmed(iSess) = sum(lc_NPH); % num NPH cells for each session
        numGi_confirmed(iSess) = sum(lc_Gi);   % num Gi cells for each session
        
        numNPH_best_guess(iSess) = sum(bg_NPH); % num NPH cells for each session
        numGi_best_guess(iSess) = sum(bg_Gi);   % num Gi cells for each session
        
        %         %% Check for Confidence Score from lesion evidence
        %
        %         if isfield(ExpKeys, 'LesionStructureConfirmed')
        %             if isfield(ExpKeys, 'LesionConfidence')        % see if 'LesionConfidence' exists in ExpKeys.  ***It should in every session.
        %                 if isempty(ExpKeys.LesionConfidence)  % This field shoule be a cell {}. See if it is empty.
        %                     lesion_conf_NPH(iSess) = NaN;  % if empty, then ExpKeys.RecrodingStructureBestGuess should be NOT empty.
        %                 else
        %                     lesion_conf_NPH(iSess) = ExpKeys.LesionConfidence;
        %                 end
        %             else
        %                 lesion_conf_NPH(iSess) = NaN; disp('ExpKeys.LesionConfidence field is not present.')
        %             end
        %         else
        %             disp('ExpKeys.LesionConfidence field is not present.')
        %         end
        %
        %
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
    else
        disp('VIDEO FILE NOT FOUND')
    end
end
Z.numNPH_confirmed = numNPH_confirmed; Z.sumNPH_confirmed = sum(Z.numNPH_confirmed);
Z.numGi_confirmed = numGi_confirmed;   Z.sumGi_confirmed = sum(Z.numGi_confirmed);

Z.numNPH_best_guess = numNPH_best_guess; Z.sumNPH_best_guess = sum(Z.numNPH_best_guess);
Z.numGi_best_guess = numGi_best_guess;   Z.sumGi_best_guess = sum(Z.numGi_best_guess);


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











