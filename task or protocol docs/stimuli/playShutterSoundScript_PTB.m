% 2022-10-28. JJS.
% Plays recorded audio (object) of the laser shutter opening and closing,
% with the same timing as that used during head-fixed opto stim.

tic;
t1 = PsychPortAudio('Start', pahandle, repetitions, 0, 1);
toc

[succeeded, cheetahReply] = NlxSendCommand('-PostEvent "ShutterSound On" 128 666');  % the variables are [Event String, TTL, Event ID] 

if succeeded == 1
    disp 'event posted'
    return;
end
