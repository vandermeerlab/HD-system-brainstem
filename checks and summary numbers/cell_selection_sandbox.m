%% set up data path
cd('C:\data\U01\datatouse');
cfg = [];
%cfg.rats = {'M085', 'M089', 'M090'}; % only specific folders
f = dir; f = f(3:end-1); f = f([f.isdir]); % grab all folder names
cfg.rats = {f.name};

fd = getDataPath(cfg);

%%
LesionConfidenceRatings = [];
for iS = 1:length(fd)

    fprintf('Entering session %d...\n', iS);
    pushdir(fd{iS});

    LoadExpKeys;

    if isempty(ExpKeys.RecordingTarget) | ~strcmp(ExpKeys.RecordingTarget{1}, 'NPH') 

        disp('Recording not in NPH, skipping...');

        popdir;
    end

    LesionConfidenceRatings = [LesionConfidenceRatings ExpKeys.LesionConfidence];

    popdir;
end