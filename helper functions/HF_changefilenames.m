function HF_changefilenames
% 2020-03-30. JJS. Code taken from RenameFiles.m.
disp('working');
fd = FindFiles('*Events.nev');
startSess = 1;
endSess = length(fd);

RenamedSomething = 0;
for iSess = startSess:endSess;
    pushdir(fileparts(fd{iSess}));
    SSNstr = HD_GetSSN('SingleSession');
    disp(SSNstr);
    %% Rename Keys File    
    Temp = dir('*keys.m');
    if exist(Temp.name, 'file')>=1 && ~strcmp(Temp.name,strcat(strrep(SSNstr, '-', '_'), '_keys.m'));
        RenamedSomething=1;
        disp('Renaming keys file');
        File=Temp.name;
        name = strcat(strrep(SSNstr, '-', '_'), '_keys.m');
        java.io.File(File).renameTo(java.io.File(name));
    end
    %% Rename events file
    Temp = dir('*.Nev');
    if exist(Temp.name, 'file')>=1 && ~strcmp(Temp.name,strcat(SSNstr, '-Events.Nev'));
        RenamedSomething =1;
        disp('Renaming events file.');
        File= Temp.name;
        java.io.File(File).renameTo(java.io.File(strcat(SSNstr, '-Events.Nev')));
    end
    
    %% Rename mp4 file 
    Temp = dir('*.mp4');
    if exist(Temp.name, 'file')>=1 && ~strcmp(Temp.name,strcat(SSNstr, '-VT1.mp4'));
        RenamedSomething =1;
        disp('Renaming mp4 file.');
        File= Temp.name;
        java.io.File(File).renameTo(java.io.File(strcat(SSNstr, '-VT1.mp4')));
    end
    
    %% Rename nvt file 
    Temp = dir('*.nvt');
    if exist(Temp.name, 'file')>=1 && ~strcmp(Temp.name,strcat(SSNstr, '-VT1.nvt'));
        RenamedSomething =1;
        disp('Renaming nvt file.');
        File= Temp.name;
        java.io.File(File).renameTo(java.io.File(strcat(SSNstr, '-VT1.nvt')));
    end    
    
    %% Rename smi file 
    Temp = dir('*.smi');
    if exist(Temp.name, 'file')>=1 && ~strcmp(Temp.name,strcat(SSNstr, '-VT1.smi'));
        RenamedSomething =1;
        disp('Renaming smi file.');
        File= Temp.name;
        java.io.File(File).renameTo(java.io.File(strcat(SSNstr, '-VT1.smi')));
    end       
    
    
%     Temp = dir('DD*.mat');
%     if exist(Temp.name, 'file')>=1 && ~strcmp(Temp.name,strcat(SSNstr, '-DD.mat'));
%         RenamedSomething =1;
%         disp('Renaming DD.mat file.');
%         File= Temp.name;
%         java.io.File(File).renameTo(java.io.File(strcat(SSNstr, '-DD.mat')));
%     end
    %% Rename TT Files
    for number=1:8
        numstr=num2str(number);
        
        %Format TT numbers one through nine with a leading zero 01, 02, ...
        if number < 10
            formatnumstr=['0' numstr];
        else
            formatnumstr=numstr;
        end
        
        
        if exist(cat(2,'Sc',numstr,'.ntt'), 'file')>=1 || exist(cat(2,'TT',numstr,'.ntt'), 'file')>=1
            fprintf(strcat('Renaming TT',formatnumstr, '\n'));
            NTTFile=[dir(cat(2,'Sc',numstr,'.ntt')) dir(cat(2,'TT',numstr,'.ntt'))];
            NTTstr=[SSNstr '-TT' formatnumstr];
            File=NTTFile(1).name;
            java.io.File(File).renameTo(java.io.File(strcat(NTTstr, '.ntt')));
        end
    end 
    %% Rename CSC Files
    %All .csc files are named Rrrrr-yyyy-mm-dd-CSCttx.ncs (where tt is 01 to 12
    %[you will need to look at the experiment sheet to find which tetrode the
    %CSC was recorded from; and x is a/b/c/d referring to the channel from
    %which the csc was recorded).
    %     Refs = {CSC1, CSC2, CSC3, CSC4, CSC5, CSC6, CSC7, CSC8, CSC9, CSC10, CSC11, CSC12, CSC13, CSC14, CSC15, CSC16 ...
    %         CSC17, CSC18, CSC19, CSC20, CSC21, CSC22, CSC23, CSC24, CSC25, CSC26, CSC27, CSC28, CSC29, CSC30, CSC31, CSC32};
    
    for number=1:36
        numstr=num2str(number);
        %Format CSC numbers one through nine with a leading zero 01, 02, ...
        if number < 10
            formatnumstr=['0' numstr];
        else
            formatnumstr=numstr;
        end
        
        %         ch = Refs{number};
        if exist(cat(2,'CSC',numstr,'.Ncs'), 'file')>=1
            RenamedSomething =1;
            fprintf(strcat('Renaming CSC',formatnumstr, '\n'));
            NCSFile=dir(cat(2,'CSC',numstr,'.Ncs'));
            NCSstr=[SSNstr '-CSC' formatnumstr];
            File=NCSFile(1).name;
            java.io.File(File).renameTo(java.io.File(strcat(NCSstr, '.Ncs')));
        end
    end
    %%
    if RenamedSomething == 0,
        disp('No files were found to be renamed.');
    end
    
    popdir;
end


