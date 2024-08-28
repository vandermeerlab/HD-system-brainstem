function [pass_flag] = CheckHFfiles(cfg_in, fd)
% 2024-07-16. JJS. For brainstem AHV/EYE neuron recording project. 
% CheckHFreqs: This function looks through each recording session in the directory and checks that the required files are present. 
%   A separate function called CheckHFkeys checks whether the critical fields are present in the keys file. 
% The manadory files (in addition to the CSC and ntt files) are ... 
%   Original files from the recording session:  *Events.nev, *VT1.mp4, *VT1.smi, *CSC21.ncs or CSC33.Ncs for platform encoder, CSC35.ncs, CSC36.Ncs for wheel encoder
%   Postprocessing eyetracking files:   *VT1_proc.mat (from Facemap), *saccades-edited.mat (semi-manual), 
%   Postprocessing tetrode files:       *.wv ClusterQual.mat  
%   Other:                              *keys.m, *AHV_StationaryTimes
%% 
if isempty(fd) 
    fd = FindFiles('*keys.m');
end
cfg_def.startSess = 1;
cfg_def.endSess = length(fd);
cfg_def.requireExpKeys = 1;  % ExpKeys is needed for every session.
cfg_def.require_mp4 = 1;     % 
cfg_def.require_smi = 1; 
cfg_def.require_proc_mat = 1;
cfg_def.require_saccades_edited = 1;
cfg_def.require_AHV_StationaryTimes = 1;
cfg_def.require_Event = 1;
cfg_def.require_ClusterQual = 1;
cfg_def.require_wv = 1;
cfg_def.require_platform_encoder = 1;
cfg_def.require_wheel = 1;

cfg = ProcessConfig(cfg_def,cfg_in);

disp([mfun,': searching for required data in session folders...'])
pass_flag = 1;

for iSess = cfg_out.startSess : cfg_out.endSess
    [path, sessID, ~] = fileparts(fd{iSess}; 
    pushdir(path)
%     SSN = HD_GetSSN; disp(SSN)
        if cfg.requireExpKeys
            fn = FindFiles('*keys.m');
            if isempty(fn)
                disp(['ExpKeys file not found in ',sessID])
                pass_flag = 0;
            end
            if ~isempty(fn) && ~isempty(cfg.ExpKeysFields)
                [ismissing,~] = checkfields(cfg,'*keys.m',cfg.ExpKeysFields);
                if ismissing
                    pass_flag = 0;
                end
            end
                    if cfg.requireVT
            fn = FindFiles('*.nvt');
            if isempty(fn)
                disp(['Video tracking file not found in ',sessionID])
                vt_fn = FindFiles('*VT1.zip');
                if ~isempty(vt_fn) & cfg.unzip
                    system(cat(2,'7z x ',vt_fn{1}));
                end
                pass_flag = 0;
            end
        end
        end
            if cfg.requireCandidates
            fn = FindFiles('*-candidates.mat');
            if isempty(fn)
                disp(['Candidates file not found in ',sessionID])
                pass_flag = 0;
            end
        end
        
        
            if ~pass_flag 
        disp(' ') % for formatting, kind of
    end
end
if pass_flag
    disp('checkTmazeReqs: all known requisites exist')
end
end




























