function [ac,xbin] = acf(spk_t,binsize,max_t)

% function [ac,xbin] = acf(spike_times,binsize,max_t)
%
% estimate autocorrelation of input spike train
%
% INPUTS:
% spike_times: [1 x nSpikes] double of spike times (in s)
% binsize: acf bin size in s
% max_t: length of acf in s
%
% OUTPUTS:
% ac: autocorrelation estimate (normalized so that acf for the zero bin is 0)
% xbin: bin centers (in s)
%
% MvdM 2013


xbin_centers = -max_t-binsize:binsize:max_t+binsize; % first and last bins are to be deleted later

ac = zeros(size(xbin_centers));

for iSpk = 1:length(spk_t)

   relative_spk_t = spk_t - spk_t(iSpk);

   ac = ac + hist(relative_spk_t,xbin_centers); % note that hist() puts all spikes outside the bin centers in the first and last bins! delete later.

end

xbin = xbin_centers(2:end-1); % remove unwanted bins

ac = ac(2:end-1);

ac = ac./max(ac); % normalize

