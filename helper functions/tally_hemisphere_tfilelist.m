function [hemisphere, left, right, numLeft, numRight, perLeft, perRight] = tally_hemisphere_tfilelist(tfilelist)

for iNeuron = 1:length(tfilelist)
    [path, neuron_to_use, ext] = fileparts(tfilelist{iNeuron});
    neuronList{iNeuron} = neuron_to_use;
    disp(neuron_to_use)
    if strcmp(pwd, path) == 0    % if current neuron is in a new folder, cd to that folder
        pushdir(path);
        SSN = HD_GetSSN;
    end
    EvalKeys
    neuronID = strcat(neuron_to_use, ext);
    cfg_spikes.fc = {neuronID};  % do one neuron at a time
    cellname{iNeuron} = neuronID; % this is the name of the eye tracking neurons (SSN-TT__.t, etc.)
    
    fn = FindFiles('*.t'); 
    [~, b, c] = fileparts(fn);
    d = strcat(b, c); 
    yesno = strcmp(cellname{iNeuron}, d);
    index = find(yesno);
    hemisphere{iNeuron} = ExpKeys.Hemisphere{index};
end

left = strcmp({'L'}, hemisphere); 
right = strcmp({'R'}, hemisphere); 
numLeft = sum(left);
numRight = sum(right);
perLeft = numLeft/length(tfilelist);
perRight = numRight/length(tfilelist);



% numLeft, numRight, perLeft, perRight