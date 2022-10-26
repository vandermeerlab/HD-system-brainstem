function ADVolts = GetADVolts(fn)

H = ReadNewHeader(fn);
L = find(~cellfun(@isempty,strfind(H, '-ADBitVolts')));
[~,remain] = strtok(H{L});
ADVolts = str2double(remain);
