function [] = HF_CreateBlankRecordingKeys(varargin)
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
    SSN = HD_GetSSN;
    disp(SSN);
    
    fout = cat(2,SSN,'_keys.m');
    fout = regexprep(fout,'-','_');
    
    does_exist = exist(fout);
    if does_exist == 2
        disp('keys file already exists. Check this file to see if you want to overwrite it.')
        return
    end
    
    fid = fopen(fout,'w');
    %populate keys.m
    fprintf(fid,'ExpKeys.Behavior = ''%s'';\n','HF Task');               % What is the task? Headfixed, freely moving, etc. 
    fprintf(fid,'ExpKeys.Protocol = ''%s'';\n','Si Probe Recording');    % What kind of recordings are being done (probe type, or type of electrode)  
    fprintf(fid,'ExpKeys.ProbeType = ''%s'';\n','A4x2tet');  
    fprintf(fid,'ExpKeys.TimeOnTrack = %d;\n',NaN);   % if there is a pre-record and post-record period, this is the first timestamp of task onset. 
    fprintf(fid,'ExpKeys.TimeOffTrack = %d;\n',NaN);   % if there is a pre-record and post-record period, this is the end of the task
    fprintf(fid,'ExpKeys.VirusTarget1 = ''%s'';\n',''); % Name of brain area for virus injection 1, if applicable 
    fprintf(fid,'ExpKeys.VirusTarget2 = ''%s'';\n',''); % Name of brain area for virus injection 2, if applicable
    fprintf(fid,'ExpKeys.RecordingTarget = ''%s'';\n',''); % Which brain regions are being recorded from. 
    fprintf(fid,'ExpKeys.TetrodeTargets = ''%s'';\n',''); % Where each tetrode is recording from (areas above).
    fprintf(fid,'ExpKeys.DepthFromSkull = %d;\n',NaN); % Estimated depth from skull surface. Based on stereotax values.  
    fprintf(fid,'ExpKeys.Notes = %d;\n', NaN);
        
    disp('keys file created')
end



