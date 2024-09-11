function [X] = cell_summary_statistics(fd)
% JJS. 2024-08-20.
% This function tabulates the number of NPH and Gi neurons from all recording sessions in the input directory, as well as location confidence scores.

% Inputs:       fd - file directory. Usually this is a list of all recording sessions with eyetracking. Ifempty, then searches in current dir.

% Outputs:

if isempty(fd)
    fd = FindFiles('*keys.m');
end
confirmed_field = NaN(1,length(fd));  % pre-allocate
best_guess_field = NaN(1,length(fd));
structure_mismatch = 0;

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
        lc_PnC = strcmp(ExpKeys.LesionStructureConfirmed, 'PnC');               % pontine reticular nucleus.            Gi turns into PNC in the anterior aspect
        lc_DMTg = strcmp(ExpKeys.LesionStructureConfirmed, 'DMTg');             % dorsomedial tegmental nucleus
        
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
        X.numPnC_confirmed(iSess) = sum(lc_PnC);
        X.numDMTg_confirmed(iSess) = sum(lc_DMTg);
        
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
        X.numPnC_confirmed(iSess) = NaN;
        X.numDMTg_confirmed(iSess) = NaN;
    end
    
    numConfirmed(iSess) = sum(X.numNPH_confirmed(iSess) + X.numGi_confirmed(iSess) + X.numSGN_confirmed(iSess) + X.numMVN_confirmed(iSess) + X.numDTN_confirmed(iSess) + ...
        X.numAbducens_confirmed(iSess) + X.num6N_confirmed(iSess) + X.nummlf_confirmed(iSess) + X.numPDTg_confirmed(iSess) + X.numLDTg_confirmed(iSess) + X.numCG_confirmed(iSess) + ...
        X.nummRt_confirmed(iSess) + X.numPnC_confirmed(iSess) + X.numDMTg_confirmed(iSess), 'omitnan');  % how many neurons present in the Confirmed field
    
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
        bg_PnC = strcmp(ExpKeys.RecordingStructureBestGuess, 'PnC');            % pontine reticular nucleus
        bg_DMTg = strcmp(ExpKeys.RecordingStructureBestGuess, 'DMTg');          % dorsomedial tegmental nucleus
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
        X.numPnC_guess(iSess) = sum(bg_PnC);
        X.numDMTg_guess(iSess) = sum(bg_DMTg);
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
        X.numPnC_guess(iSess) = NaN;
        X.numDMTg_guess(iSess) = NaN;
        X.numunknown_guess(iSess) = NaN;
    end
    
    numGuess(iSess) = sum(X.numNPH_guess(iSess) + X.numGi_guess(iSess) + X.numSGN_guess(iSess) + X.numMVN_guess(iSess) + X.numDTN_guess(iSess) + ...
        X.numAbducens_guess(iSess) + X.num6N_guess(iSess) + X.nummlf_guess(iSess) + X.numPDTg_guess(iSess) + X.numLDTg_guess(iSess) + X.numCG_guess(iSess) + ...
        X.nummRt_guess(iSess) + X.numPnC_guess(iSess) + X.numDMTg_guess(iSess) + X.numunknown_guess(iSess), 'omitnan');  % how many neurons present in the BestGuess field
    
    if numGuess(iSess) > 0; guess_structureID(iSess) = 1; end   % 0 indicates no confirmed neurons. 1 indicates some confirmed neuron(s).
    
    if confirmed_structureID + guess_structureID ~= 1
        warning('problem with structureID keys fields');  % One of these fields should have a strucutre ID for one or more neurons. Both fields should NOT have entries in the same keys file.
    end
    
    tfiles = FindFiles('*.t');
    numCells = length(tfiles);
    %     assert(numConfirmed(iSess) + numGuess(iSess) == numCells);
    if numConfirmed(iSess) + numGuess(iSess) ~= numCells                     % Check to make sure that each neuron has a structureID entry in keys file.
        warning('# neurons does not match # of structureID entries')
        structure_mismatch = structure_mismatch + 1;
    end
    
    %% Check for Confidence Score from lesion evidence
    if isfield(ExpKeys, 'LesionConfidence')        % see if 'LesionConfidence' exists in ExpKeys.  ***It should in every session.
        lesion_confidence_field(iSess) = 1;
        if isempty(ExpKeys.LesionConfidence)  % This field shoule be a cell {}. See if it is empty.
            lesion_confidence_score(iSess) = NaN;  % if empty, then ExpKeys.RecrodingStructureBestGuess should be NOT empty.
        else
            lesion_confidence_score(iSess) = ExpKeys.LesionConfidence;
        end
    else
        lesion_confidence_field(iSess) = NaN;
        lesion_confidence_score(iSess) = Inf; 
        warning('ExpKeys.LesionConfidence field is not present.')
    end
    %% Check for Confidence Score from best guess evidence
    if isfield(ExpKeys, 'NeuronConfidence')        % see if 'LesionConfidence' exists in ExpKeys.  ***It should in every session.
        neuron_confidence_field(iSess) = 1;
        if isempty(ExpKeys.NeuronConfidence)  % This field shoule be a cell {}. See if it is empty.
            neuron_confidence_score(iSess) = NaN;  % if empty, then ExpKeys.RecrodingStructureBestGuess should be NOT empty.
        else
            neuron_confidence_score(iSess) = ExpKeys.NeuronConfidence;
        end
    else
        neuron_confidence_field(iSess) = NaN;
        neuron_confidence_score(iSess) = Inf; 
        warning('ExpKeys.NeuronConfidence field is not present.')
    end
    %%
    clear ExpKeys
    popdir;
    clear confirmed_structureID
    clear guess_structureID
end

%% Calculate the total for each structure for each session 
X.NPH = X.numNPH_confirmed + X.numNPH_guess; 
X.Gi = X.numGi_confirmed + X.numGi_guess; 
X.SGN = X.numSGN_confirmed + X.numSGN_guess; 
X.MVN = X.numMVN_confirmed + X.numMVN_guess; 
X.DTN = X.numDTN_confirmed + X.numDTN_guess; 
X.Abducens = X.numAbducens_confirmed + X.numAbducens_guess; 
X.s6N = X.num6N_confirmed + X.num6N_guess; 
X.mlf = X.nummlf_confirmed + X.nummlf_guess; 
X.PDTg = X.numPDTg_confirmed + X.numPDTg_guess; 
X.LDTg = X.numLDTg_confirmed + X.numLDTg_guess; 
X.CG = X.numCG_confirmed + X.numCG_guess; 
X.mRt = X.nummRt_confirmed + X.nummRt_guess; 
X.PnC = X.numPnC_confirmed + X.numPnC_guess; 
X.unknown = X.numunknown_guess; 

%% Tally how many recording sessions apply to each brain area 
X.num_sess_NPH = length(find(X.NPH));
X.num_sess_Gi = length(find(X.Gi));
X.num_sess_SGN = length(find(X.SGN));
X.num_sess_MVN = length(find(X.MVN));
X.num_sess_DTN = length(find(X.DTN));
X.num_sess_Abducens = length(find(X.Abducens));
X.num_sess_6N = length(find(X.s6N));
X.num_sess_mlf = length(find(X.mlf));
X.num_sess_PDTg = length(find(X.PDTg));
X.num_sess_LDTg = length(find(X.LDTg));
X.num_sess_CG = length(find(X.CG));
X.num_sess_mRt = length(find(X.mRt));
X.num_sess_PnC = length(find(X.PnC));
X.num_sess_unknown = length(find(X.unknown));

%% Calculate the sum total over sessions -- CONFIRMED
X.sumNPH_confirmed = sum(X.numNPH_confirmed, 'omitnan');
X.sumGi_confirmed = sum(X.numGi_confirmed, 'omitnan');
X.sumSGN_confirmed = sum(X.numSGN_confirmed, 'omitnan');
X.sumMVN_confirmed = sum(X.numMVN_confirmed, 'omitnan');
X.sumDTN_confirmed = sum(X.numDTN_confirmed, 'omitnan');
X.sumAbducens_confirmed = sum(X.numAbducens_confirmed, 'omitnan');
X.sum6N_confirmed = sum(X.num6N_confirmed, 'omitnan');
X.summlf_confirmed = sum(X.nummlf_confirmed, 'omitnan');
X.sumPDTg_confirmed = sum(X.numPDTg_confirmed, 'omitnan');
X.sumLDTg_confirmed = sum(X.numLDTg_confirmed, 'omitnan');
X.sumCG_confirmed = sum(X.numCG_confirmed, 'omitnan');
X.summRt_confirmed = sum(X.nummRt_confirmed, 'omitnan');
X.sumPnC_confirmed = sum(X.numPnC_confirmed, 'omitnan');

%% Calculate the sum total over sessions -- BEST GUESS
X.sumNPH_guess = sum(X.numNPH_guess, 'omitnan');
X.sumGi_guess = sum(X.numGi_guess, 'omitnan');
X.sumSGN_guess = sum(X.numSGN_guess, 'omitnan');
X.sumMVN_guess = sum(X.numMVN_guess, 'omitnan');
X.sumDTN_guess = sum(X.numDTN_guess, 'omitnan');
X.sumAbducens_guess = sum(X.numAbducens_guess, 'omitnan');
X.sum6N_guess = sum(X.num6N_guess, 'omitnan');
X.summlf_guess = sum(X.nummlf_guess, 'omitnan');
X.sumPDTg_guess = sum(X.numPDTg_guess, 'omitnan');
X.sumLDTg_guess = sum(X.numLDTg_guess, 'omitnan');
X.sumCG_guess = sum(X.numCG_guess, 'omitnan');
X.summRt_guess = sum(X.nummRt_guess, 'omitnan');
X.sumPnC_guess = sum(X.numPnC_guess, 'omitnan');
X.sumunknown_guess = sum(X.numunknown_guess, 'omitnan');

%% Calculate the total over both categories -- CONFIRMED + BESTGUESS
X.NPH_total = X.sumNPH_confirmed + X.sumNPH_guess;
X.Gi_total = X.sumGi_confirmed + X.sumGi_guess;
X.SGN_total = X.sumSGN_confirmed + X.sumSGN_guess;
X.MVN_total = X.sumMVN_confirmed + X.sumMVN_guess;
X.DTN_total = X.sumDTN_confirmed + X.sumDTN_guess;
X.Abducens_total = X.sumAbducens_confirmed + X.sumAbducens_guess + X.sum6N_confirmed + X.sum6N_guess;
X.mlf_total = X.summlf_confirmed + X.summlf_guess;
X.PDTg_total = X.sumPDTg_confirmed + X.sumPDTg_guess;
X.LDTg_total = X.sumLDTg_confirmed + X.sumLDTg_guess;
X.CG_total = X.sumCG_confirmed + X.sumCG_guess;
X.mRt_total = X.summRt_confirmed + X.summRt_guess;
X.PnC_total = X.sumPnC_confirmed + X.sumPnC_guess;
X.unknown_total = X.sumunknown_guess;

X.total_neurons1 = X.NPH_total + X.Gi_total + X.SGN_total + X.MVN_total + X.DTN_total + X.Abducens_total + ...
    X.mlf_total + X.PDTg_total + X.LDTg_total + X.CG_total + X.mRt_total + X.PnC_total + X.unknown_total;

X.num_confirmed = numConfirmed;
X.total_num_confirmed = sum(numConfirmed);
X.num_guess = numGuess;
X.total_num_guess = sum(numGuess);
X.num_all = X.num_confirmed + X.num_guess;
X.structure_mismatch = structure_mismatch;

X.total_neurons2 = X.total_num_confirmed + X.total_num_guess;

if X.total_neurons1 ~= X.total_neurons2; warning('neuron count mismatch'); end

X.confirmed_field = confirmed_field;    % this is brain structure ID for each (confirmed) neuron 
X.guess_field = best_guess_field;       % this is brain structure ID for each (best guess) neuron 

X.lesion_confidence_field = lesion_confidence_field; % NaN = field does not exist. 0 = field is empty. 1 = field has a value.  This is the confidence score for confirmed neurons. 
X.neuron_confidence_field = neuron_confidence_field; 

% Assert that neuron confidence entries and lesion confidence entries are complementary sets and do not overalp
X.lesion_confidence_score = lesion_confidence_score;   % confidence scores for lesion confirmed session.              NaN = no entry. 1 = high confidence. 2 = moderate confidence. 3 = uncertain. 
X.neuron_confidence_score = neuron_confidence_score;   % confidence scores for best guess sessions                    NaN = no entry. 1 = high confidence. 2 = moderate confidence. 3 = uncertain. 

X.confirmed_sessions = find(~isnan(X.lesion_confidence_score));
X.best_guess_sessions = find(~isnan(X.neuron_confidence_score));

X.problem_sessions = intersect(X.confirmed_sessions, X.best_guess_sessions); 
if ~isempty(X.problem_sessions); warning('overlap detected between confidence scores for Confirmed and BestGuess sessions'); end 

if length(X.confirmed_sessions) + length(X.best_guess_sessions) ~= length(fd); warning('mismatch between number sessions in dir and the number of sessions with confidence field'); end
%% Print the main results









% temp_confirmed = X.sumNPH_confirmed + X.sumGi_confirmed + X.sumSGN_confirmed + X.sumMVN_confirmed + X.sumDTN_confirmed + X.sumAbducens_confirmed + ...
%     X.sum6N_confirmed + X.summlf_confirmed + X.sumPDTg_confirmed + X.sumLDTg_confirmed + X.sumCG_confirmed + X.summRt_confirmed + X.sumPnC_confirmed; % 189
%
% temp_guess = X.sumNPH_guess + X.sumGi_guess + X.sumSGN_guess + X.sumMVN_guess + X.sumDTN_guess + X.sumAbducens_guess + ...
%     X.sum6N_guess + X.summlf_guess + X.sumPDTg_guess + X.sumLDTg_guess + X.sumCG_guess + X.summRt_guess + X.sumPnC_guess + ...  % 51
%     X.sumunknown_guess;

