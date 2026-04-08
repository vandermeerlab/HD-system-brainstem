function [tt_place, sessCounter, simultaneous_NPH_cells] = find_tt_nums2(tfilelist)
% JJS. 2024-10-01. This function searches through the neurons in tfilelist to determine how many cells were recorded from each tetrode (over all cells)

% Input:        tfilelist - a cell array of .t files (with their path)

% Output:       tt_place - a list of tetrode numbers for all the recorded neurons. Should be same length as number of neurons in tfilelist.
doPlot = 0;
sessCounter = 0;
cellcounter = 0;
% simultaneous_NPH_cells = NaN(1,length(tfilelist));
if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
for iNeuron = 1:length(tfilelist)
    [path, neuron_to_use, ~] = fileparts(tfilelist{iNeuron});
    if strcmp(pwd, path) == 0
        pushdir(path);
        EvalKeys
        disp(neuron_to_use)
        sessCounter = sessCounter + 1;
        tS = FindFiles('*.t');
        [~, neuronIDs, ~] = fileparts(tS);
        if ~iscell(neuronIDs); neuronIDs = {neuronIDs}; end
        index = strcmp(neuron_to_use, neuronIDs);  % this result should be a nCell x 1 logical. if [0 1 0] for instance, then the second neuron in S is our spike train to use.
        row_to_use = find(index);
        
        lc_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'NPH');
        bg_NPH = strcmp(ExpKeys.RecordingStructureBestGuess, 'NPH');            % search for NPH string for best guess recordings
        if isempty(lc_NPH) % || sum(lc_NPH) == 0
            temp1 = zeros(1,length(tS))';
            lc_NPH = logical(temp1);
        end
        if isempty(bg_NPH) % || sum(bg_NPH) == 0
            temp2 = zeros(1,length(tS))';
            bg_NPH = logical(temp2);
        end
        comb_NPH = lc_NPH + bg_NPH; % combined array of zeros and ones to indicate if each neuron is NPH.
        simultaneous_NPH_cells(sessCounter) = length(find(comb_NPH));
    end
    ind = strfind(neuronIDs, 'TT');
    for iCelltoUse = 1:length(row_to_use)
        cellcounter = cellcounter + 1;
        tt_place(cellcounter) = str2num(neuronIDs{row_to_use(iCelltoUse)}(ind{iCelltoUse}+3));
    end
end

if doPlot
    clf
    edges = [1:8];
    hist(tt_place,edges)
    xlabel('tetrode #')
    ylabel('# of neurons')
    set(gca, 'FontSize', 18)
    c = axis;
    axis([c(1) 9 c(3) c(4)])
end

popdir;
disp(num2str(sum(simultaneous_NPH_cells))); 


