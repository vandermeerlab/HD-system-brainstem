function [nasal_indices_REST, temporal_indices_REST, nasal_timestamps_REST, temporal_timestamps_REST, numNasal, numTemporal, nNaNtouse, tNaNtouse, sessList] = isolateStationarySaccadesAllSessions(fd, varargin)

% JJS. 2022-09-05.
% This function finds the indices (which) and timestamps (when) of saccades that fall during stationary periods when the platform is at rest.
% Uses the function isolateStationarySaccades.m to count the saccades within a single session folder. This function gathers data for all sessions, and organizes it.

% Inputs:       fd -     list of sessions to run through. Can be empty, in which case it will run through the current directory.
% Outputs:      nasal_indices  -    1 x nSession cell array where each cell is one session. For nasal saccades only.
%                                       integers -      indicates the indices of which saccades from the matlab var 'nasalSaccades' are from the stationary period.
%                                       empty cells -   indicates sessions with stationary periods, but no saccades that occured during those epochs.
%                                       NaNs -          indicates sessions in which no video data exists (camera not set up yet)
%
%               temporal_indices -  same as above, for temporal saccades

%               nasalREST     -     1 x nSession cell array. Same as above. these values indicate saccade timestamps. For nasal saccades.
%               temporalREST  -     1 x nSession cell array. Same as above. these values indicate saccade timestamps. For temporal saccades.

%               numNasal      -     1 x nSession double.
%                                       zero -          indicates that zero saccasdes occured during stationary periods.
%                                       NaN  -          indicates that no video data exists for this session.
%                                       integers -      indicates how many saccades occured during stationary periods.
%               numTemporal   -     1 x nSession double. Same as above.
%               nNaNtouse     -     binary array which indicates which sessions had video tracking (1 = yes. 0 = no).
%               tNaNtouse     -     binary array which indicates which sessions had video tracking (1 = yes. 0 = no). For temporal saccades
%               doPlot        -     option to plot the saccades (evoked and stationary).

doPlot = 1;
MarkerSize = 10;
FontSize = 18;
saccadeThresh = 15; 
process_varargin(varargin); 

if isempty(fd) 
    fd = FindFiles('*keys.m');
end
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; disp(SSN);
    %     if exist(strcat(SSN, '-VT1_proc.mat')) ==2
    %         disp('video file found')
%     [nasal_indices{iSess}, temporal_indices{iSess}, nasalREST{iSess}, temporalREST{iSess}] = isolateStationarySaccades;
    [nasal_indices_REST{iSess}, temporal_indices_REST{iSess}, nasal_timestamps_REST{iSess}, temporal_timestamps_REST{iSess}] = isolateStationarySaccades();

    %     else
    %         disp('video file not detected. Skipping.')
end
% popdir;
% end
%% replace NaNs with empty brackets in order to determine the length of each element.
% ind = cellfun(@isnan, nasal_indices, 'UniformOutput', false);

nasal_indices_noNaN = nasal_indices_REST;
temporal_indices_noNaN = temporal_indices_REST;
for iSess = 1:length(nasal_indices_REST)
    if isnan(nasal_indices_REST{iSess})
        nasal_indices_noNaN{iSess} = [];
    end
    if isnan(temporal_indices_REST{iSess})
        temporal_indices_noNaN{iSess} = [];
    end
    % keep track of sessions without video data [versus session with video data, but no stationary saccades]
    nNaN{iSess} = isnan(nasal_indices_REST{iSess});
    if find(nNaN{iSess})
        nNaNtouse(iSess) = 0;   % this indicates a session with ___NO___ video tracking data
    else
        nNaNtouse(iSess) = 1;   % this indicates a session with ___YES___ video tracking data
    end
    tNaN{iSess} = isnan(temporal_indices_REST{iSess});
    if find(tNaN{iSess})
        tNaNtouse(iSess) = 0;   % this indicates a session with ___NO___ video tracking data
    else
        tNaNtouse(iSess) = 1;   % this indicates a session with ___YES___ video tracking data
    end
end
% re-insert NaNs
numNasal = cellfun(@length, nasal_indices_noNaN);  % NaNs were removed prior to this step becasue length([NaN]) = 1;
numNasal(nNaNtouse==0) = NaN; % re-insert the NaNs for sessions with no video data.
numTemporal = cellfun(@length, temporal_indices_noNaN);
numTemporal(tNaNtouse==0) = NaN; % re-insert the NaNs for sessions with no video data.

whichSess = numNasal >= saccadeThresh & numTemporal >= saccadeThresh;
sessList = fd(whichSess);
sessList = fileparts(sessList); 
numSess = sum(whichSess);
disp('Number of sessions with at least the minimum number of saccades')
disp(num2str(numSess)) 


if doPlot == 1
    plot(numNasal, 'MarkerSize', MarkerSize); hold on
    plot(numTemporal, 'MarkerSize', MarkerSize)
    xlabel('Session Number', 'FontSize', FontSize)
    ylabel('Number of Stationary Saccades', 'FontSize', FontSize)
    title('Spontaneous Saccades', 'FontSize', FontSize)
    set(gca, 'XTick', [1:length(fd)])
end



end