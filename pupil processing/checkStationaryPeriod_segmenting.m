function checkStationaryPeriod_segmenting(fd, numNasal, numTemporal)


for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; 
    if isnan(numNasal(iSess)) || isnan(numTemporal(iSess))
        disp('skipping session. no video data')
    else
%         if exist(strcat(SSN, '-saccades-edited.fig')) == 2
            open(strcat(SSN, '-saccades-edited.fig'))
            plotSTART_STOP_lines;
            disp(strcat('session number ', num2str(iSess)))
            disp(strcat('num nasal stationary = ', num2str(numNasal(iSess))))
            disp(strcat('num temporal stationary = ', num2str(numTemporal(iSess))))
            pause
%             clf
%         end
    end
    popdir;
end

    