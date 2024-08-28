function [pass_flag] = CheckHFkeys(cfg_in)
% 2024-08-28. JJS. For brainstem AHV/EYE neuron recording project. 
% CheckHFreqs: This function looks through each recording session in the directory and checks that the required files and ExpKey fields are present. 
%   Detailed explanation goes here

% cfg_def.requireExpKeys = 0;
% cfg_def.ExpKeysFields = {}; list of strings specifying field names
%           ex: {'nTrials','badTrials'}
% cfg_def.requireMetadata = 0;
% cfg_def.MetadataFields = {}; list of strings specifying field names
% cfg_def.requireVT = 0;
% cfg_def.requireCandidates = 0;
% cfg_def.requireTimes = 0; for R042 only
% cfg_def.requireHSdetach = 0; for R044 only
% cfg_def.requireFiles = 0; files folder for images generated during
%   analysis
% cfg_def.ratsToProcess = {'R042','R044','R050','R064'}; % only process these
%   rats
% cfg_def.verbose = 1; 1 display command window text, 0 don't
% cfg_def.unzip = 1; % attempt to unzip vt file if no *.nvt found with 7zip.org

%% 
cfg_def.checkall = 0;

cfg_def.requireExpKeys = 1;  % ExpKeys is needed for every session.
cfg_def.ExpKeysFields = {};  % List of all of the essential keys fields. This list is not exhaustive. The ExpKeys fields have changed over time (usually by addition). 
cfg_def.require_mp4 = 1;     % 
cfg_def.require_mp4 = 1;
cfg_def.require_smi = 1; 
cfg_def.require_proc_mat = 1;
cfg_def.require_saccades_edited = 1;
cfg_def.require_AHV_StationaryTimes = 1;
cfg_def.require_Event = 1;
cfg_def.require_ClusterQual = 1;
cfg_def.require_wv = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);
cfg.verbose = 1; % i want checkfields to tell me which ones are missing

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

