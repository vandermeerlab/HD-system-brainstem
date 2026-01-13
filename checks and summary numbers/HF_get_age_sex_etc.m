function [numSess, sessID, sex, age] = HF_get_age_sex_etc()

cd('C:\Jeff\U01\datatouse')
load('NPH all tfilelist.mat') %#ok<LOAD>
numSess = length(NPH_all_tfilelist); %#ok<USENS>

for iSess = 1:length(NPH_all_tfilelist)
    [a, b, c] = fileparts(NPH_all_tfilelist{iSess});
    pushdir(fileparts(NPH_all_tfilelist{iSess}));
    SSN = HD_GetSSN; disp(SSN)
    EvalKeys;
    sessID{iSess} = b;
    sex{iSess} = ExpKeys.sex;
%     if ~isempty(ExpKeys.age)
        age(iSess) = ExpKeys.age;
%     else
%         age(iSess) = NaN;
%     end
    popdir;
end

age = age';
sex = sex';
sessID = sessID';