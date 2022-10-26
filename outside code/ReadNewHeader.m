function H = ReadNewHeader(fn)

fp = fopen(fn, 'r');
if fp==-1
    error('Cannot open %s.\n', fn);
end
H = {};
while ftell(fp) < 16384
    H{end+1} = fgetl(fp);
end
end

