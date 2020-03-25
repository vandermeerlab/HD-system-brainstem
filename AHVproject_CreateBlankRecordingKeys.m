function [] = AHVproject_CreateBlankRecordingKeys(varargin)
%% This function Generates the Keys File with the optional inputs left blank

fd = FindFiles('*Events.nev');
startSess = 1;
endSess = length(fd);

process_varargin(varargin);
for iSess = startSess:endSess;
    pushdir(fileparts(fd{iSess}));
    
    %%%%%%%%%%%%%%%%%%%
    % generate keys.m %
    %%%%%%%%%%%%%%%%%%%
    SSN = PM_GetSSN;
    disp(SSN);
    
    fout = cat(2,SSN,'_keys.m');
    fout = regexprep(fout,'-','_');
    fid = fopen(fout,'w');
    %populate keys.m
    fprintf(fid,'ExpKeys.Behavior = ''%s'';\n','PM Task');
    fprintf(fid,'ExpKeys.Protocol = ''%s'';\n','Behavior');
    
    fprintf(fid,'ExpKeys.NorthArmFlavor = ''%s'';\n','Banana');
    fprintf(fid,'ExpKeys.WestArmFlavor = ''%s'';\n','Chocolate');
    fprintf(fid,'ExpKeys.SouthArmFlavor = ''%s'';\n','Banana');
    fprintf(fid,'ExpKeys.EastArmFlavor = ''%s'';\n','Chocolate');
    
    fprintf(fid,'ExpKeys.DevalFlavor = ''%s'';\n',''); % Options are 'Grain' or 'Banana'
    
    fprintf(fid,'ExpKeys.numPelletsNorth = %d;\n',3);
    fprintf(fid,'ExpKeys.numPelletsWest = %d;\n',3);
    fprintf(fid,'ExpKeys.numPelletsSouth = %d;\n',1);
    fprintf(fid,'ExpKeys.numPelletsEast = %d;\n',1);
    
    fprintf(fid,'ExpKeys.NorthPelletsLeftOver = %d;\n',NaN);
    fprintf(fid,'ExpKeys.WestPelletsLeftOver = %d;\n',NaN);
    fprintf(fid,'ExpKeys.SouthPelletsLeftOver = %d;\n',NaN);
    fprintf(fid,'ExpKeys.EastPelletsLeftOver = %d;\n',NaN);
    
    fprintf(fid,'ExpKeys.Weight = %0.1f;\n',NaN);
    fprintf(fid,'ExpKeys.PostFeed = %0.1f;\n', NaN);
    
    fprintf(fid,'ExpKeys.TimeOnTrack = %d;\n',NaN);
    fprintf(fid,'ExpKeys.TimeOffTrack = %d;\n',NaN);
    fprintf(fid,'ExpKeys.TaskPhase = ''%s'';\n',''); % Options here include acquistion, overtraining, probe session, post-deval
    fprintf(fid,'ExpKeys.TrainingGroup = ''%s'';\n',''); %  Options include 'CT' (for trained to criterion) or 'OT' (overtrained)
    fprintf(fid,'ExpKeys.VirusGroup = ''%s'';\n',''); %
    fprintf(fid,'ExpKeys.PostDevalCons = ''%s'';\n','');
    fprintf(fid,'ExpKeys.Notes = ''%s'';\n','');
    
    fprintf(fid,'ExpKeys.Target = %d;\n', NaN);
    fprintf(fid,'ExpKeys.TetrodeTargets = %d;\n', NaN');
    
    disp('keys file created')
end



% fprintf(fid,'ExpKeys.Behavior = ''%s'';\n',Behavior);
% fprintf(fid,'ExpKeys.Protocol = ''%s'';\n',Protocol);
% fprintf(fid,'ExpKeys.PelletRatio = %d;\n',PelletRatio);
% fprintf(fid,'ExpKeys.Tones = %d;\n',Tones);
% fprintf(fid,'ExpKeys.Weight = %0.1f;\n',Weight);
% fprintf(fid,'ExpKeys.Blocks = %d;\n', Blocks);
% fprintf(fid,'ExpKeys.Nudges = %d;\n', Nudges);
% fprintf(fid,'ExpKeys.PostFeed = %0.1f;\n', PostFeed);
% fprintf(fid,'ExpKeys.FeederR1 = %d;\n', RF);
% fprintf(fid,'ExpKeys.FeederL1 = %d;\n', LF);
% fprintf(fid,'ExpKeys.TimeOnTrack = %d;\n',TimeOnTrack);
% fprintf(fid,'ExpKeys.TimeOffTrack = %d;\n',TimeOffTrack);
% fprintf(fid,'ExpKeys.Notes = ''%s'';\n',Notes);
