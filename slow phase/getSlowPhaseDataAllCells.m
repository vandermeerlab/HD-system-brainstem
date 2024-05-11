function [] = getSlowPhaseDataAllCells(cfg_in)

% JJS. 2024-04-17. This function runs getSlowPhaseData on all sessions in the directory and saves the plots from each neuron into the specified location.


fd = FindFiles('*keys.m');
cfg_def = [];
cfg_def.doSave = 1;
cfg_def.OneLocation = true;
cfg_def.startSess = 1;
cfg_def.endSess = length(fd);
cfg_master = ProcessConfig2(cfg_def, cfg_in);

if cfg_master.OneLocation
    fpath = 'D:\Jeff\U01\analysis\dot mat files\SlowPhasePlots'; % save all figs to on location
else
    fpath = pwd; % save each figure to its original session folder
end

for iSess = cfg_master.startSess : cfg_master.endSess
    disp(iSess)
    pushdir(fileparts(fd{iSess}));
    % Skip over sessions that don't have eye tracking data
    SSN = HD_GetSSN;
    if exist(strcat(SSN, '-VT1.smi')) == 2
        sd = LoadSessionData([]);
        for iCell = 1:length(sd.S.t)
            figure('units','normalized','outerposition',[0 0 1 1])
            [~, ~] = getSlowPhaseData(cfg_in, sd, iCell);
            if cfg_master.doSave == 1
                f = fullfile(fpath, strcat(sd.fn{iCell}, '-slowphasedata.jpg'));
                saveas(gcf, f)
                disp('data saved')
            end
        end
        popdir;
    else
        disp('no video file detected. Skipping session.')
        popdir;
    end
end