function ts_out = ts_timeOffset(ts_in, time_offset)
% function ts_out = ts_timeOffset(ts_in, time_offset)
%
% add time_offset to every cell in ts_in
%
% common example: if you got the start time of a tsd from somewhere
% like [~, start_time] = tsd_startAtZero() or csc.tvec(1) and you want to apply 
% that time to a ts so that they have the same start time, then you should 
% do (note minus sign):
%
% ts_new = ts_timeOffset(ts_old, -start_time)

if ~CheckTS
    error('Input not a correctly formed ts.');
end
   
add_offset_fun = @(x) x + time_offset;
ts_out = ts_in;
ts_out.t = cellfun(add_offset_fun, ts_out.t, 'UniformOutput', false);