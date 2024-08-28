function [pass_flag] = CheckHFkeys2(cfg_in)
% 2024-08-28. JJS. For brainstem AHV/EYE neuron recording project. 
% CheckHFreqs: This function looks through each recording session in the directory and checks that the required files and ExpKey fields are present. 
%   The exact version of the keys file has changed many times over the course of the headfixed brainstem project. The latter recording sessions tend to
%   have more keys fields. There are some fields that do need to be in all keys files and those are emphasized here. I will add the option in this code
%   to add those elements that are not present to existing keys files and also to check for other, optional fields that are not necessarily critical. 

% List of all of the fields to check for. Most important ones are listed first. 

% ver2. This version has a simpler structure and lists all the required keys fields in a single cell array. 
cfg_def.require_ExpKeys.notes = 0; % {}; General notes about the session. Text could be anything here.  
cfg_def.require_ExpKeys.eventLabels = 0;  % {}; list of strings that indicates the identity of each event code. Could include 'optical stimulation', 'lights out recording', etc. 
cfg_def.require_ExpKeys.eventNumbers = 0; % []; list of numerical event IDs in the Events.Nev file that correspond to the ExpKeys.eventLabels. This will be same length as ExpKeys.eventLabels. 
cfg_def.require_ExpKeys.eventMeaning = 0; % {};        % for instance, {'Laser On', 'Laser Off'} 
cfg_def.require_ExpKeys.experimenter = 0; % {'jstott'};   % person who conducted the recording itself 
cfg_def.require_ExpKeys.species = 0; % {}; 'Mouse';      % 'Mouse' or 'Rat' 
cfg_def.require_ExpKeys.sex =  0; % [];           % 'Male' or 'Female'
cfg_def.require_ExpKeys.strain =  0; % [];    % name of the mouse strain. 'C57BL/6J' (WT) or 'Ai32' or 'CHETA'.   
cfg_def.require_ExpKeys.JAXnum =  0; % [];            % strain number, if using a JAX mosue. 000664 for wild-type. 24109 for Ai-32.  
cfg_def.require_ExpKeys.project =  0; % 'AHV brainstem';    % AHV brainstem = recording in brainstem (could be NPH/SGN/DTN/DPGi/MVN). 'AHV ADN' = ADN 
cfg_def.require_ExpKeys.phase =  0; % 'performance';
cfg_def.require_ExpKeys.behavior =  0; % 'passive rotation';
cfg_def.require_ExpKeys.EncoderCSC =  0; % []; 
cfg_def.require_ExpKeys.Quadrature1CSC =  0; % [];
cfg_def.require_ExpKeys.Quadrature2CSC =  0; % [];
cfg_def.require_ExpKeys.AnalogEncoder =  0; % []; 
cfg_def.require_ExpKeys.VirusTarget1 =  0; % {'DTN'};  % if applicable. For mouse with no virus injection, leave blank. 
cfg_def.require_ExpKeys.VirusTarget2 =  0; % {};
cfg_def.require_ExpKeys.VirusID_1 =  0; % {'pENN.AAV.hSyn.HI.eGFP-Cre.WPRE.SV40'}; % 'pENN.AAV.hSyn.HI.eGFP-Cre.WPRE.SV40' for retroAAV. 
cfg_def.require_ExpKeys.VirusID_1stocknum =  0; % 105540;  % Addgene plasmid number, if applicable
cfg_def.require_ExpKeys.VirusID_2 =  0; % {};   % example   {'pAAV-EF1a-double floxed-hChR2(H134R)-mCherry-WPRE-HGHpA'}
cfg_def.require_ExpKeys.VirusID_2stocknum =  0; % []; % Addgene plasmid number, if applicable
cfg_def.require_ExpKeys.RecordingDay =  0; % []; % what day of recording i.e. how many times have I recorded previously, in terms of actually recorded sessions).
cfg_def.require_ExpKeys.DaysPostSurgery =  0; % []; % days since virus injection surgery 
cfg_def.require_ExpKeys.DaysPostCraniotomy =  0; % []; % days since the craniotomy was performed
cfg_def.require_ExpKeys.weight =  0; % []; % weight in grams
cfg_def.require_ExpKeys.bday =  0; % []; % mouse date of birth  % YYYY-MM-DD. Not a string. 
cfg_def.require_ExpKeys.age =  0; % []; % in days 
cfg_def.require_ExpKeys.ProbeNum =  0; % '';  % ID number for the probe. 
cfg_def.require_ExpKeys.Plating =  0; % {}; % is there plating of the probe? If not, leave empty. Common strings will be 'IrOx' (iridium oxide) and 'Z-coat'. 
cfg_def.require_ExpKeys.ProbeMaker =  0; % {'NN'};   % Can be 'CNT' for Cambridge Neurotech or 'NN' for Neuronexus. 
cfg_def.require_ExpKeys.ProbeType =  0; % {'A4x2-tet'}; 
cfg_def.require_ExpKeys.ProbeChannels =  0; % 1:32;
cfg_def.require_ExpKeys.RecordingTarget =  0; % {'NPH'};
cfg_def.require_ExpKeys.RecordingTarget2 =  0; % {'Left', 'Right'};
cfg_def.require_ExpKeys.ProbeFiber =  0; % [];   % Was an optical fiber present on the probe? value of 0 or 1. 0 = 'no'. 1 = 'yes'. NaN = unknown
cfg_def.require_ExpKeys.LaserOn =  0; []; % was the laser used in the recording; i.e., was opto stimulation attempted (apart from whether it was succesful)?
cfg_def.require_ExpKeys.FiberDiameter =  0; % [];  % in microns
cfg_def.require_ExpKeys.ShankWidth =  0; % []; % how wide are the probe shanks? 15 microns or 50 microns are typical values. 
cfg_def.require_ExpKeys.FiberNA =  0; %  [];  % numerica aperature. For example, 0.22 NA
cfg_def.require_ExpKeys.FiberPower =  0; % []; % estimated light output, in milliwatts (mW). If unknown, leave blank. 
cfg_def.require_ExpKeys.IgnoreOpto =  0; % []; % If a value of 1 is added here, that means that even if this session has a cell(s) that pass the threshold for being opto responsive, they should be ignored. A value of 1 would be appropriate if histology indicates that the virus missed DTN or that the recording might be in some area other than DTN projecting (for whatever reason).  
cfg_def.require_ExpKeys.OpticalResponse =  % 0; []; % Was there an obvious, subjective optical response? Answer 0 for 'no' and 1 for 'yes'. 
cfg_def.require_ExpKeys.DIOlabel =  0; % 0; % logical value of 0 or 1 to indicate whether DIO was painted onto the probe before recording. 
cfg_def.require_ExpKeys.DIOcolor =  0; % {}; % what color is the dye (i.e. 'red' or 'green' or 'far_red')
cfg_def.require_ExpKeys.DIOpresent =  0; % []; % was the DIO label present in the histology for this session? can add explanation in ExpKeys.Notes. 0 = 'no'. 1 = 'yes'; Leave blank or use NaN if unknown for whatever reason.
cfg_def.require_ExpKeys.MarkingLesion.made =  0; % []; % logical value of 0 or 1 to indicate whether a current was passed through the probe to make a marking lesion. Use a value of 2 to indicate that a lesion from a different recording, but in the same electrode penetration is visible. 
cfg_def.require_ExpKeys.Lesion.present =  0; []; % is a lesion even present which could correspond to this session? 0 = 'no'. 1 = 'yes'. If 'yes', then fill out ExpKeys.LesionConfidence below.  Use a value of 2 to indicate that a lesion from a different recording, but in the same electrode penetration is visible. Yes here could indicate that shanks are visible, even if lesion mark is not. If so, indicate in the notes.
cfg_def.require_ExpKeys.MarkingLession.channel =  0; % []; % which channels was/were used. ADchannel. Empty for no marking lesion made.
cfg_def.require_ExpKeys.MarkingLession.tetrode =  0; % [];  % which tetrode was used for the marking lesion??? Empty for no marking lesion made.
cfg_def.require_ExpKeys.MarkingLession.tetrodechannel =  0; % []; % which channel of the tetrode was used? 1-4
cfg_def.require_ExpKeys.MarkingLession.current =  0; % []; % how much current was used (in microamps). Empty for no marking lesion made.
cfg_def.require_ExpKeys.MarkingLession.duration =  0; % []; % how long was current applied (in seconds). Empty for no marking lesion made.
cfg_def.require_ExpKeys.MarkingLession.polarity =  0; % {'unipolar'}; % typically this is unipolar. The options are 'unipolar' or 'bipolar'.
cfg_def.require_ExpKeys.MarkingLession.directionality =  0; % {'cathodal'}; % options here are 'cathodal' or 'anodal'. Cathodal = red lead attached to probe -- recording electrode site. black lead should be attached to animal ground. The red light on the stimulus isolator should be over the red wire (and tape that says 'cathodal')
cfg_def.require_ExpKeys.MarkingLession.reps =  0; % []; % How many times was the stimulus repeated? If only applied once, enter '1' or leave blank. For multiple reps, enter value.  
cfg_def.require_ExpKeys.LesionSiteNotes =  0; % {};
cfg_def.require_ExpKeys.EyeCamera =  0; % 1; % was a camera present for eye tracking? 0 = 'no'. 1 = 'yes'.
cfg_def.require_ExpKeys.EyeTrackingQuality =  0; % []; % value from 1-5 indicating quality of eye tracking, if an eyetracking file is present. If not present, leave blank. '5' indicates best. '3' is ok. '1' is not usable.
cfg_def.require_ExpKeys.EyeTrackingNotes =  0; % {}; 
%% Critical Section relating to histology and electrode location and placement 
% LOCATION LOCATION LOCATION
cfg_def.require_ExpKeys.AnteriorPosition =  0; % []; % bregma value for the anterior/posterior position of the recording plane, as estimated from histology, or measurements from skull landmark. If no iformation, then leave empty.   
cfg_def.require_ExpKeys.MedialLateral =  0; % [NaN]; % estimated and/or intended M/L position of the center of the recording probe. This will almost always be zero. 
cfg_def.require_ExpKeys.DorsalVentral =  0; % [NaN]; % estimated and/or intended D/V position of the center of the recording probe. This would come from histology, if different from ExpKeys.ProbeDepth 
cfg_def.require_ExpKeys.probeDepth =  0; % []; % this will likely be the D/V number as read from the stereotax at the time of recording. For verified lesion sites, it should be the number indicated in the 'matlas' Excel spreadsheet. 
cfg_def.require_ExpKeys.LocationNotes =  0; % {};
cfg_def.require_ExpKeys.FluorescenceNotes =  0; % {};
cfg_def.require_ExpKeys.ThioninNotes =  0; % {}; 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% FOR NEURONS where we have LESIONS
cfg_def.require_ExpKeys.LesionStructureConfirmed =  0; % {}; % This should be a string which designates which structure was recorded from, if it is obvious from the histology.
% Appropriate entries include 'NPH', 'SGN', 'DTN', 'Abducens', 'Gi', 'mlf', 'DMTg', 'MVN', 'PDTg, 'LDTg', 'CG' (central grey), 'PnC' (pontine reticular nucleus). For multiple neurons, enter strings into the cell, separated by semicolons. 
cfg_def.require_ExpKeys.LesionConfidence =  0; % [];  % based on histology (if it exists), how confident are we about where the lesion mark is/which lesion mark is which? 1 = Extremely confident. 2 = somewhat confident. 3 = unsure. 
%Empty for no marking lesion made. 0 = current was applied but no mark could be found. This will be a single numeral 
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% FOR NEURONS where we are estimating position without lesions 
cfg_def.require_ExpKeys.RecordingStructureBestGuess =  0; % {}; % This is the best guess of which structure was recorded from, if it can be determined from notes, approximate location, etc. 
    % this entry can remain empty if "ExpKeys.LesionStructureConfirmed" is filled out. like {'6N'; '6N'; 'Gi'; 'mlf'; 'PNC' (pontine reticular nucleus)}; For multiple neurons, enter strings into the cell, separated by semicolons. 
    % Appropriate entries include '6N', 'NPH', 'SGN', 'DTN', 'Abducens', 'Gi', 'mlf', 'DMTg', 'MVN', 'PDTg, 'LDTg', 'CG' (central grey), 'mRt' (mesenphalic reticular formation). For multiple neurons, enter strings into the cell, separated by semicolons. 
%ExpKeys.NeuronCoordinates = {}; % matrix with the probe's estimated medial/lateral (M/L) and dorsal/ventral (D/V) coordiantes from rough estimates when histological LESIONS are not present. ***This field may be unnecessary, as it is included in the matlas file. 
cfg_def.require_ExpKeys.NeuronConfidence =  0; % []; % (mostly, the A/P position) namely, how confident in the brain region and rough estimate of position. 1 = confident. 2 = some confidence. 3 = low confidence. If no information, then leave empty. If ExpKeys.NeuronCoordinates is empty, then NeuronConfidence will be empty. 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Which HEMISPHERE
cfg_def.require_ExpKeys.Hemishpere =  0; % {}; % this array should be the same length as the number of neurons, and should specify 'L' (left), 'R' (right) or 'C' (center) for each neuron. 'Center' applies to mlf area where the lesion is likely dead center. If totally unsure, leave the whole array empty 
% or include a 'NaN' value for a particular neuron that is unclear. This should be filled out whether the location method is "lesion" or "educated guess".
cfg_def.require_ExpKeys.HemisphereConfidence =  0; % [];   % Based on observations during recording and histology, how confident are we about which hemisphere each TT was in? 1 = very confident. 2 = somewhat confident. 3 = not confident. Empty = no information. 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cfg = ProcessConfig(cfg_def,cfg_in,mfun);


% cfg_def.ExpKeysFields = {};  % List of all of the essential keys fields. This list is not exhaustive. The ExpKeys fields have changed over time (usually by addition). 

if cfg.checkall
    cfg.requireExpKeys = 1;
    cfg.ExpKeysFields = {'RestrictionType','Session','Layout','Pedestal','pathlength','patharms','realTrackDims','convFact','nPellets','waterVolume','nTrials','forcedTrials','nonConsumptionTrials','badTrials','TimeOnTrack','TimeOffTrack','prerecord','task','postrecord','goodSWR','goodTheta'};
    cfg.requireMetadata = 1;
    cfg.MetadataFields = {'coord','taskvars','SWRtimes','SWRfreqs'};
    cfg.requireVT = 1;
    cfg.requireCandidates = 1;
    cfg.requirePrecandidates = 1;
    cfg.requireTimes = 1; % R042 only
    cfg.requireHSdetach = 1; % R044 only
    cfg.requireFiles = 1;
end
%% 
disp([mfun,': searching for required data in session folders...'])
filesep = '\';

%get data path
base_fp = getBaseFP;

%%

pass_flag = 1;

% Remember where you started
original_folder = pwd; % pwd is your current directory ("print working directory")

end

%% Old keys field names
% ExpKeys.IrOx = 1;  %  0 = 'no' and 1 = 'yes'. This is iridium oxide plating, offered by Neuronexus. Increases charge density and reduces risk of probe damage with stim. If unknown, leave blank.  


