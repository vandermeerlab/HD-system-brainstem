fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}))
    SSN = HD_GetSSN; disp(SSN);
    HF_changefilenames;
    popdir;
end

threshold = 2; 
for iCell = 1: size(zN, 1);
    ZN(iCell,:) = zN(iCell,:)
    
    zN(zN > 2) = 1;  
    
    Xn = zN; 
    Xn(zN > threshold) = 1;
    Xn(zN < -threshold) = -1; 
    loggy = zN < threshold & zN > -threshold; 
    Xn(loggy) = 0; 
    
    Xt = zT; 
    Xt(zT > threshold) = 1;
    Xt(zT < -threshold) = -1; 
    loggy = zT < threshold & zT > -threshold; 
    Xt(loggy) = 0; 
    
    
end
    %%
    
    fd = FindFiles('*keys.m');
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}))
    SSN = HD_GetSSN; disp(SSN);
    if exist(strcat(SSN, '-VT1.smi')) == 2
        assert(exist(strcat(SSN, '-VT1_proc.mat')) == 2)
    end
    popdir;
end


for iCell = 1:length(FRxBinN);
    dN(iCell,:) = diff(FRxBinN(iCell,:));
    dT(iCell,:) = diff(FRxBinT(iCell,:));
end
