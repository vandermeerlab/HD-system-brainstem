function [] = 

fd = FindFiles('*.mp4');
cfg_in.plot_time_series = 0;
cfg_in.plot_PETH = 0;
CONTROL_out = [];
PUFF_out = [];
for iSess = 1:length(fd)
    pushdir(fileparts(fd{iSess}));
    [puff_out, control_out, cfg_out] = airpuff_peth(cfg_in);
    
    PUFF_out = vertcat(PUFF_out, puff_out.data);
    CONTROL_out = vertcat(CONTROL_out, control_out.data);
end
