function convert_temporal_nasal_to_ipsi_contra(fd,cfg_in)
%2024-09-20. JJS. This function looks at every .keys file and adds IPSI & CONTRA labels for each neuron. 
%   (alongside the existing tmeporal and nasal). All saccade variables (e.g., timestamps & amplitudes) are therefore saved with both temporal-nasal and 
%   ispi-contra labels. The data itelf is identical. The ipsi - contra distinction is determined by which brain hemisphere the neuron is in (w/ respect to saccade direction).

%   Key to variable naming. ALWAYS FOR THE LEFT EYE! 

%    LEFT hemisphere neuron, temporal (leftward) saccade   = IPSI 
%    LEFT hemisphere neuron, nasal (rightward) saccade      = CONTRA 

%    RIGHT hemisphere neuron, temporal (leftward) saccade  = CONTRA 
%    RIGHT hemisphere neuron, nasal (rightward) saccade     = IPSI 


cfg_out = ProcessConfig(cfg_def, cfg_in);

if isempty(fd)
    fd = FindFiles('*keys.m');
end
numSess = length(fd);
for iCell = 1:numSess
    pushdir(fileparts(fd{iSess}));
    SSN = HD_GetSSN; disp(SSN)
    load(strcat(SSN, '-saccades-edited.mat'))
    EvalKeys;








end

