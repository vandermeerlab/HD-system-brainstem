function [tc_wheel] = getWheel_TC(sd, cfg_in)
% JJS. 2023-04-28.
% Calculate and plot the tuning curve for eye position
% input:   sd - session data structure with spike trains S
% output:  pupilTC - a 1 x nCell array of structures, each with the fields tc, occ_hist, spk_hist, usr, cfg 

