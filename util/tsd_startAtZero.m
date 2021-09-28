function [tsd_out, first_t] = tsd_startAtZero(tsd_in)
% function [tsd_out, first_t] = tsd_startAtZero(tsd_in)
%
% resets a tsd's tvec to start at zero by subtracting tsd.tvec(1), returns
% the result as tsd_out and the first timestamp as first_t

if ~CheckTSD(tsd_in)
   error('Input is not a tsd'); 
end

first_t = tsd_in.tvec(1);

if first_t == 0
   warning('First timestamp is already zero. That''s an unlikely coincidence.'); 
end

tsd_out = tsd_in;
tsd_out.tvec = tsd_out.tvec - first_t;