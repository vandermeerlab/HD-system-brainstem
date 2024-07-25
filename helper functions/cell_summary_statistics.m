function [numNPH, numGi, conf_NPH, Z, numGuessNPH, NPHguess_conf] = cell_summary_statistics(fd)

if isempty(fd)
    fd = FindFiles('*keys.m');
end
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}))
    SSN = HD_GetSSN; disp(SSN)
    if exist(strcat(SSN,'-VT1.mp4'))
        EvalKeys;
        tf_NPH = strcmp(ExpKeys.LesionStructureConfirmed, 'NPH');
        tf_Gi = strcmp(ExpKeys.LesionStructureConfirmed, 'Gi');
        numNPH(iSess) = sum(tf_NPH);
        numGi(iSess) = sum(tf_Gi);
        if isfield(ExpKeys, 'LesionConfidence')
            if isempty(ExpKeys.LesionConfidence)
                conf_NPH(iSess) = NaN;
            else
                conf_NPH(iSess) = ExpKeys.LesionConfidence;
            end
        else
            conf_NPH(iSess) = NaN;
        end
        if isfield(ExpKeys, 'RecordingStructureBestGuess')
            if isempty(ExpKeys.RecordingStructureBestGuess)
                numGuessNPH(iSess) = NaN;
            else
                guess_NPH = strcmp(ExpKeys.RecordingStructureBestGuess, 'NPH');
                numGuessNPH(iSess) = sum(guess_NPH);
                if isempty(ExpKeys.NeuronConfidence)
                    warning('this ExpKeys should have a neuron confidence value. However, it is blank')
                    NPHguess_conf(iSess) = NaN;
                else
                NPHguess_conf(iSess) = ExpKeys.NeuronConfidence;
                end
            end
        else
            NPHguess_conf(iSess) = NaN;
        end
    else
        numNPH(iSess) = NaN;
        numGi(iSess) = NaN;
        numGuessNPH(iSess) = NaN;
    end
end

Z.a = conf_NPH ==1;
Z.b = numNPH(Z.a);
Z.c = sum(Z.b);

Z.d = conf_NPH ==2;
Z.e = numNPH(Z.d);
Z.f = sum(Z.e);

Z.h = conf_NPH ==3;
Z.i = numNPH(Z.h);
Z.j = sum(Z.i);


% ExpKeys.RecordingStructureBestGuess = {};
