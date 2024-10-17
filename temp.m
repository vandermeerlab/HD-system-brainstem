function [] = temp(tfilelist)

% JJS. 2024-10-02. Calculates the saccade peth for each neuron in the tfilelist and generates a combined heatplot of saccade response for all neurons.
%                  This function is in temporal-nasal space. A separate function uses IPSI-CONTRA space. 

doPlot = 0;
sessCounter = 0;
cellcounter = 0;
% simultaneous_NPH_cells = NaN(1,length(tfilelist));
if isempty(tfilelist)
    tfilelist = FindFiles('*.t');
end
for iNeuron = 1:length(tfilelist)
    [path, neuron_to_use, ~] = fileparts(tfilelist{iNeuron});
    if strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder 
        pushdir(path);
        SSN = HD_GetSSN; disp(SSN)
        EvalKeys
        disp(neuron_to_use)
        sessCounter = sessCounter + 1;
        tS = FindFiles('*.t');
        [~, neuronIDs, ~] = fileparts(tS);
        
    end
    cfg_peth.FRwindow = [-.2 .2];
    cfg_peth.dt = 0.005;
%     load(strcat(SSN, '-saccades-edited.mat'))
    Z{iNeuron} = makeSaccadeHeatPlot(cfg_peth, [], []);