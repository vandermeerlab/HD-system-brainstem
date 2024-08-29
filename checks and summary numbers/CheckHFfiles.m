function [X] = CheckHFfiles(cfg_in, fd)
% 2024-07-16. JJS. For brainstem AHV/EYE neuron recording project.
% CheckHFreqs: This function looks through each recording session in the directory and checks that the required files are present.
%   A separate function called CheckHFkeys checks whether the critical fields are present in the keys file.
% The manadory files (in addition to the CSC and ntt files) are ...
%   Original files from the recording session:  *Events.nev, *VT1.mp4, *VT1.smi, *CSC21.ncs or CSC33.Ncs for platform encoder, CSC35.ncs, CSC36.Ncs for wheel encoder
%   Postprocessing eyetracking files:   *VT1_proc.mat (from Facemap), *saccades-edited.mat (semi-manual),
%   Postprocessing tetrode files:       *.wv ClusterQual.mat
%   Other:                              *keys.m, *AHV_StationaryTimes
%
%   Inputs:         cfg_in  = optional config settings. If you want to check all, then set cfg_in.checkall to 1.
%   Outputs:        X       = structure with fields that enumerate for each file type which sessions are missing. X.pass_flag = NaN for all files found and 0 for a session with any missing file type.
if isempty(fd)
    fd = FindFiles('*keys.m');
end
cfg_def.startSess = 1;
cfg_def.endSess = length(fd);
cfg_def.requireExpKeys = 1;  % ExpKeys is needed for every session.
cfg_def.requireExpKeysFields = 0; % require specific keys fields, as defined later in cfg.ExpKeysFields
cfg_def.require_mp4 = 1;     %
cfg_def.require_smi = 1;
cfg_def.require_proc_mat = 1;
cfg_def.require_saccades_edited = 1;
cfg_def.require_AHV_StationaryTimes = 1;
cfg_def.require_Events = 1;
cfg_def.require_ClusterQual = 1;
cfg_def.require_wv = 1;
cfg_def.require_platform_encoder = 1;
cfg_def.require_wheel = 1;
cfg_def.ExpKeysFields = {'probeDepth','LesionStructureConfirmed','LesionConfidence','RecordingStructureBestGuess','NeuronConfidence','Hemishpere',...
    'HemisphereConfidence','MarkingLesion.made','MarkingLesion.present','MarkingLession.channel','MarkingLession.tetrode','MarkingLession.current',...
    'MarkingLession.duration','MarkingLession.polarity','MarkingLession.directionality','MarkingLession.reps','MarkingLLesionSiteNotes',...
    'AnteriorPosition'};

cfg = ProcessConfig(cfg_def,cfg_in);

if ~isempty(cfg_in) && cfg_in.checkall
    cfg.requireExpKeys = 1;
    cfg.ExpKeysFields = {'probeDepth','LesionStructureConfirmed','LesionConfidence','RecordingStructureBestGuess','NeuronConfidence','Hemishpere',...
        'HemisphereConfidence','MarkingLesion.made','MarkingLesion.present','MarkingLession.channel','MarkingLession.tetrode','MarkingLession.current',...
        'MarkingLession.duration','MarkingLession.polarity','MarkingLession.directionality','MarkingLession.reps','MarkingLLesionSiteNotes',...
        'AnteriorPosition'};
    cfg.requireExpKeys = 1;  % ExpKeys is needed for every session.
    cfg.require_mp4 = 1;     %
    cfg.require_smi = 1;
    cfg.require_proc_mat = 1;
    cfg.require_saccades_edited = 1;
    cfg.require_AHV_StationaryTimes = 1;
    cfg.require_Events = 1;
    cfg.require_ClusterQual = 1;
    cfg.require_wv = 1;
    cfg.require_platform_encoder = 1;
    cfg.require_wheel = 1;
end
cfg.startSess = 1;
cfg.endSess = length(fd);
X.pass_flag = NaN(1,length(fd));

keyscounter1 = 1; keyscounter2 = 1; vidcount = 1; smi_count = 1; proc_count = 1; stationary_count = 1; cq_count = 1; wv_count = 1; platform_count = 1; ...
    wheel_count = 1;
for iSess = cfg.startSess : cfg.endSess
    if ~isempty(fd)
        [path, ~, ~] = fileparts(fd{iSess});
        pushdir(path);
        SSN = HD_GetSSN;
    end
%     disp([SSN,': searching for required data in session folders...'])
    tfiles = FindFiles('*.t'); num_t = length(tfiles);
    %% Check keys
    if cfg.requireExpKeys
        fn = FindFiles('*keys.m');
        if isempty(fn)
            disp(['ExpKeys file not found: ',SSN])
            X.keys_absent{keyscounter1} = SSN;
            X.pass_flag(iSess) = 0;
        end
        if cfg_def.requireExpKeysFields
            if ~isempty(fn) && ~isempty(cfg.ExpKeysFields)
                [ismissing,~] = checkfields(cfg,'*keys.m',cfg.ExpKeysFields);
                if ismissing
                    disp(['keys file is incomplete: ',SSN])
                    X.keystochange{keyscounter2} = SSN;
                    X.pass_flag(iSess) = 0;
                    keyscounter2 = keyscounter2 + 1;
                end
            end
        end
    end
    %% Check Video
    if cfg.require_mp4
        fn = FindFiles('*.mp4');
        if isempty(fn)
            disp(['Eye tracking video file not found: ',SSN])
            X.missing_video{vidcount} = SSN;
            X.pass_flag(iSess) = 0;
            vidcount = vidcount + 1;
        end
    end
    %% Check SMI (video timestamps)
    if cfg.require_smi
        fn = FindFiles('*.smi');
        if isempty(fn)
            disp(['smi timestamps file not found: ',SSN])
            X.missing_smi{smi_count} = SSN;
            X.pass_flag(iSess) = 0;
            smi_count = smi_count + 1;
        end
    end
    %% Check for proc_mat Facemap output file
    if cfg.require_proc_mat
        fn = FindFiles('*proc.mat');
        if isempty(fn)
            disp(['Facemap _proc.mat file not found: ',SSN])
            X.missing_proc{proc_count} = SSN;
            X.pass_flag(iSess) = 0;
            proc_count = proc_count + 1;
        end
    end
    %% Check for Stationary times file (manual list of rotation start/stop times)
    if cfg.require_AHV_StationaryTimes
        fn = FindFiles('*AHV_StationaryTimes.mat');
        if isempty(fn)
            disp(['Stationary times file not found: ',SSN])
            X.missing_video{stationary_count} = SSN;
            X.pass_flag(iSess) = 0;
            stationary_count = stationary_count + 1;
        end
    end
    %% Check for Neuralynx Events file
    if cfg.require_Events
        fn = FindFiles('*Events.nev');
        if isempty(fn)
            disp(['Events file not found: ',SSN])
            X.missing_events{events_count} = SSN;
            X.pass_flag(iSess) = 0;
            proc_count = proc_count + 1;
        end
    end
    %% Check for tetrode Cluster Quality files
    if cfg.require_ClusterQual
        fn = FindFiles('*ClusterQual.mat');
        if isempty(fn) || length(fn) ~= num_t
            disp(['ClustQual file(s) not found or incorrect number: ',SSN])
            X.missing_cqfile{cq_count} = SSN;
            X.pass_flag(iSess) = 0;
            cq_count = cq_count + 1;
        end
    end
    %% Check for tetrode wave files
    if cfg.require_wv
        fn = FindFiles('*wv.mat');
        if isempty(fn) || length(fn) ~= num_t
            disp(['Wave file(s) not found or incorrect number: ',SSN])
            X.missing_wv{wv_count} = SSN;
            X.pass_flag(iSess) = 0;
            wv_count = wv_count + 1;
        end
    end
    %% Check for Platform encoder file
    if cfg.require_platform_encoder
        fn21 = FindFiles('*CSC21.Ncs'); fn33 = FindFiles('*CSC33.Ncs'); % platform encoder is tied in directly to either CSC21 (early sessions) or CSC33. *find switch date
        if isempty(fn21) && isempty(fn33)
            disp(['Platform encoder CSC file(s) not found: ',SSN])
            X.missing_platform_file{platform_count} = SSN;
            X.pass_flag(iSess) = 0;
            platform_count = platform_count + 1;
        end
    end
    %% Check for Wheel data
    if cfg.require_wheel
        fn35 = FindFiles('*CSC35.Ncs'); fn36 = FindFiles('*CSC36.Ncs'); % platform encoder is tied in directly to either CSC21 (early sessions) or CSC33. *find switch date
        if isempty(fn35) && ~isempty(fn36)
            disp(['Wheel encoder CSC file(s) not found: ',SSN])
            X.missing_wheel{wheel_count} = SSN;
            X.pass_flag(iSess) = 0;
            wheel_count = wheel_count + 1;
        end
    end
%     if X.pass_flag(iSess)
%         disp('checkHF files: all known requisites exist')
%     end
%     else
%         disp('one or more elements are missing')
%     end
end
