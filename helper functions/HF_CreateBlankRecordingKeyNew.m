function [] = HF_CreateBlankRecordingKeyNew(varargin)
%% This function Generates the Keys File with the optional inputs left blank
%% JJS. 2023-08-21. Updated to reflect the new standard format of the keys file. Many new fields added or altered.

fd = FindFiles('*Events.nev');
startSess = 1;
endSess = length(fd);

process_varargin(varargin);
for iSess = startSess:endSess;
    pushdir(fileparts(fd{iSess}));
    
    %%%%%%%%%%%%%%%%%%%
    % generate keys.m %
    %%%%%%%%%%%%%%%%%%%
    SSN = HD_GetSSN;
    disp(SSN);
    
    fout = cat(2,SSN,'_keys.m');
    fout = regexprep(fout,'-','_');
    
    does_exist = exist(fout);
    if does_exist == 2
        while(1)
            disp('keys file already exists. Check this file to see if you want to overwrite it.')
            m = input('Do you want to continue, y/n [y]:','s');
            if m == 'n'
                disp('stopping')
                return
            else
                break;
            end
        end
    end
    disp('some fields are populated with default values. make sure to double check')
    fid = fopen(fout,'w');
    %populate keys.m
    %     fprintf(fid,'ExpKeys.Behavior = ''%s'';\n','HF Task');               % What is the task? Headfixed, freely moving, etc.
    %     fprintf(fid,'ExpKeys.Protocol = ''%s'';\n','Si Probe Recording');    % What kind of recordings are being done (probe type, or type of electrode)
    %     fprintf(fid,'ExpKeys.ProbeType = ''%s'';\n','A4x2tet');
    %     fprintf(fid,'ExpKeys.TimeOnTrack = %d;\n',NaN);   % if there is a pre-record and post-record period, this is the first timestamp of task onset.
    %     fprintf(fid,'ExpKeys.TimeOffTrack = %d;\n',NaN);   % if there is a pre-record and post-record period, this is the end of the task
    %     fprintf(fid,'ExpKeys.VirusTarget1 = ''%s'';\n',''); % Name of brain area for virus injection 1, if applicable
    %     fprintf(fid,'ExpKeys.VirusTarget2 = ''%s'';\n',''); % Name of brain area for virus injection 2, if applicable
    %     fprintf(fid,'ExpKeys.RecordingTarget = ''%s'';\n',''); % Which brain regions are being recorded from.
    %     fprintf(fid,'ExpKeys.TetrodeTargets = ''%s'';\n',''); % Where each tetrode is recording from (areas above).
    %     fprintf(fid,'ExpKeys.DepthFromSkull = %d;\n',NaN); % Estimated depth from skull surface. Based on stereotax values.
    %     fprintf(fid,'ExpKeys.Notes = %d;\n', NaN);
    %%
    fprintf(fid,'ExpKeys.notes = ''%s'';\n','');

    fprintf(fid,'eventLabels = ''%s'';\n','{}'); % list of strings that indicates the identity of each event code. Could include 'optical stimulation', 'lights out recording', etc.
    %     ExpKeys.eventLabels = {};         % list of strings that indicates the identity of each event code. Could include 'optical stimulation', 'lights out recording', etc.
    %
    %
    %     ExpKeys.eventNumbers = [];        % list of numerical event IDs in the Events.Nev file that correspond to the ExpKeys.eventLabels.
    %
    %
    %     ExpKeys.experimenter = {'jstott'};   % person who conducted the recording itself
    %
    %
    %     ExpKeys.species = 'Mouse';      % 'Mouse' or 'Rat'
    %
    %
    %     ExpKeys.sex = [];           % 'Male' or 'Female'
    %
    %
    %     ExpKeys.strain = [];    % name of the mouse strain. 'C57BL/6J' (WT) or 'Ai32' or 'CHETA'.
    %
    %
    %     ExpKeys.JAXnum = [];            % strain number, if using a JAX mosue. 000664 for wild-type. 24109 for Ai-32.
    %
    %
    %     ExpKeys.project = 'AHV brainstem';    % AHV brainstem = recording in brainstem (could be NPH/SGN/DTN/DPGi/MVN)
    %
    %
    %     ExpKeys.phase = 'performance';
    %
    %
    %     ExpKeys.behavior = 'passive rotation';
    %
    %
    %     ExpKeys.EncoderCSC = [];
    %
    %
    %     ExpKeys.Quadrature1CSC = [];
    %
    %
    %     ExpKeys.Quadrature2CSC = [];
    %
    %
    %     ExpKeys.AnalogEncoder = [];
    %
    %
    %     ExpKeys.VirusTarget1 = {'DTN'};  % if applicable. For mouse with no virus injection, leave blank.
    %
    %
    %     ExpKeys.VirusTarget2 = {};
    %
    %
    %     ExpKeys.VirusID_1 = {'pENN.AAV.hSyn.HI.eGFP-Cre.WPRE.SV40'}; % 'pENN.AAV.hSyn.HI.eGFP-Cre.WPRE.SV40' for retroAAV.
    %
    %
    %     ExpKeys.VirusID_1stocknum = 105540;  % Addgene plasmid number, if applicable
    %
    %
    %     ExpKeys.VirusID_2 = {};   % example   {'pAAV-EF1a-double floxed-hChR2(H134R)-mCherry-WPRE-HGHpA'}
    %
    %
    %     ExpKeys.VirusID_2stocknum = 20297; % Addgene plasmid number, if applicable
    %
    %
    %     ExpKeys.RecordingDay = []; % what day of recording i.e. how many times have I recorded previously, in terms of actually recorded sessions).
    %
    %
    %     ExpKeys.DaysPostSurgery = []; % days since headbar surgery/craniotomy
    %
    %
    %     ExpKeys.weight = []; % weight in grams
    %
    %
    %     ExpKeys.bday = []; % mouse date of birth  % YYYY-MM-DD. Not a string.
    %
    %
    %     ExpKeys.age = []; % in days
    %
    %
    %     ExpKeys.ProbeNum = [];  % ID number for the probe.
    %
    %
    %     ExpKeys.IrOx = 1;  %  0 = 'no' and 1 = 'yes'. This is iridium oxide plating, offered by Neuronexus. Increases charge density and reduces risk of probe damage with stim.
    %
    %
    %     ExpKeys.ProbeMaker = {'NN'};   % Can be 'CNT' for Cambridge Neurotech or 'NN' for Neuronexus.
    %
    %
    %     ExpKeys.ProbeType = {'A4x2-tet'};
    %
    %
    %     ExpKeys.ProbeChannels = 1:32;
    %
    %
    %     ExpKeys.RecordingTarget = {'NPH'};
    %
    %
    %     ExpKeys.RecordingTarget2 = {'Left', 'Right'};
    %
    %
    %     ExpKeys.probeDepth = []; % in mm ventral from SKULL surface. This field replaces ExpKeys.tetrodeDepths. This number is the number on the stereotax arm during recording, NOT from histology.
    %
    %
    %     ExpKeys.ProbeFiber = 1;   % Was an optical fiber present on the probe? value of 0 or 1. 0 = 'no'. 1 = 'yes'.
    %
    %
    %     ExpKeys.LaserOn = []; % was the laser used in the recording; i.e., was opto stimulation attempted (apart from whether it was succesful)?
    %
    %
    %     ExpKeys.FiberDiameter = [];  % in microns
    %
    %
    %     ExpKeys.FiberNA = .22;  % numerica aperature. For example, 0.22 NA
    %
    %
    %     ExpKeys.FiberPower = []; % estimated light output, in milliwatts (mW)
    %
    %
    %     ExpKeys.DIOlabel = 0; % logical value of 0 or 1 to indicate whether DIO was painted onto the probe before recording.
    %
    %
    %     ExpKeys.DIOcolor = []; % what color is the dye (i.e. 'red' or 'green' or 'far_red')
    %
    %
    %     ExpKeys.DIOpresent = []; % was the DIO label present in the histology for this session? can add explanation in ExpKeys.Notes. 0 = 'no'. 1 = 'yes';
    %
    %
    %     ExpKeys.MarkingLesion.made = 1; % logical value of 0 or 1 to indicate whether a current was passed through the probe to make a marking lesion
    %
    %
    %     ExpKeys.Lesion.present = 1; % is a lesion even present which could correspond to this session? 0 = 'no'. 1 = 'yes'. If 'yes', then fill out ExpKeys.LesionConfidence below.
    %
    %
    %     ExpKeys.MarkingLession.channel = []; % which channels was/were used. ADchannel. Empty for no marking lesion made.
    %
    %
    %     ExpKeys.MarkingLession.tetrode = [];  % which tetrode was used for the marking lesion??? Empty for no marking lesion made.
    %
    %
    %     ExpKeys.MarkingLession.current = []; % how much current was used (in microamps). Empty for no marking lesion made.
    %
    %
    %     ExpKeys.MarkingLession.duration = []; % how long was current applied (in seconds). Empty for no marking lesion made.
    %
    %
    %     ExpKeys.MarkingLession.reps = []; % How many times was the stimulus repeated? If only applied once, enter '1' or leave blank. For multiple reps, enter value.
    %
    %     %% EYE Tracking
    %     ExpKeys.EyeCamera = 1; % was a camera present for eye tracking? 0 = 'no'. 1 = 'yes'.
    %
    %     ExpKeys.EyeTrackingQuality = []; % value from 1-5 indicating quality of eye tracking, if an eyetracking file is present. If not present, leave blank. '5' indicates best. '3' is ok. '1' is not usable.
    %
    %
    %     ExpKeys.EyeTrackingNotes = {};
    %
    %     %% LOCATION LOCATION LOCATION
    %     ExpKeys.AnteriorPosition = []; % bregma value for the anterior/posterior position of the recording plane, as estimated from histology, or measurements from skull landmark. If no iformation, then leave empty.
    %
    %
    %     ExpKeys.MedialLateral = 0; % estimated and/or intended M/L position of the center of the recording probe. This will almost always be zero.
    %
    %
    %     ExpKeys.probeDepth = []; % this will likely be the D/V number as read from the stereotax at the time of recording. For verified lesion sites, it should be the number indicated in the 'matlas' Excel spreadsheet.
    %
    %
    %     %% FOR NEURONS where we have LESIONS
    %     ExpKeys.LesionStructureConfirmed = {}; % size = 1 x n neurons. Each entry should be a semicolon-separated string which designates which structure was recorded from, if it is obvious from the histology.
    %     % Appropriate entries include '6N', 'NPH', 'SGN', 'DTN', 'Abducens', 'Gi', 'mlf', 'DMTg', 'MVN', 'PDTg, 'LDTg', 'CG' (central grey). For multiple neurons, enter strings into the cell, separated by semicolons.
    %
    %
    %     ExpKeys.LesionConfidence = [];  % based on histology (if it exists), how confident are we about where the lesion mark is/which lesion mark is which? 1 = Extremely confident. 2 = somewhat confident. 3 = unsure.
    %     %Empty for no marking lesion made. 0 = current was applied but no mark could be found. This will be a single numeral
    %
    %     %% FOR NEURONS where we are estimating position without lesions
    %     ExpKeys.RecordingStructureBestGuess = {}; % This is the best guess of which structure was recorded from, if it can be determined from notes, approximate location, etc.
    %     % this entry can remain empty if "ExpKeys.LesionStructureConfirmed" is filled out. like {'6N'; '6N'; 'Gi'; 'mlf'; 'PNC' (pontine reticular nucleus)}; For multiple neurons, enter strings into the cell, separated by semicolons.
    %     %ExpKeys.NeuronCoordinates = {}; % matrix with the probe's estimated medial/lateral (M/L) and dorsal/ventral (D/V) coordiantes from rough estimates when histological LESIONS are not present
    %
    %
    %     ExpKeys.NeuronConfidence = []; % vector with confidence values for the location of the recording; namely, how confident in the brain region and rough estimate of position. 1 = confident. 2 = some confidence. 3 = low confidence. If no information, then leave empty. If ExpKeys.NeuronCoordinates is empty, then NeuronConfidence will be empty.
    %
    %
    %     %% Which HEMISPHERE
    %     ExpKeys.Hemishpere = {}; % this array should be the same length as the number of neurons, and should specify 'L' (left), 'R' (right) or 'C' (center) for each neuron. 'Center' applies to mlf area where the lesion is likely dead center. If totally unsure, leave the whole array empty
    %     % or include a 'NaN' value for a particular neuron that is unclear. This should be filled out whether the location method is "lesion" or "educated guess".
    %
    %
    %     ExpKeys.HemisphereConfidence = [];   % One numeral for each neuron. Based on observations during recording and histology, how confident are we about which hemisphere each TT was in? 1 = very confident. 2 = somewhat confident. 3 = not confident. Empty = no information.
    %
    %
    %     ExpKeys.ThioninNotes = {};
    %
    
    
    
    
    
    disp('keys file created')
end



