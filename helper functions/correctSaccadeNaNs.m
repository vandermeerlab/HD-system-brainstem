SSN = HD_GetSSN; disp(SSN);
load(FindFile('*saccades-edited.mat'))

keep = find(~isnan(temporalSaccades));
temporalSaccades = temporalSaccades(keep);  %#ok<*FNDSB>
temporalAmplitudes = temporalAmplitudes(keep);
numtemporalSaccades = length(temporalSaccades);

keep = find(~isnan(nasalSaccades));
nasalSaccades = nasalSaccades(keep);
nasalAmplitudes = nasalAmplitudes(keep);
numnasalSaccades = length(nasalSaccades);

combinedSaccades = numtemporalSaccades + numnasalSaccades;

save(strcat(SSN, '-saccades-edited.mat'))
disp('file saved')
