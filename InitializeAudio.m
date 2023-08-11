cd('C:\Users\admin\Documents\Sound recordings')
disp(pwd)

connected = NlxAreWeConnected();
if connected == 1
    disp 'Cheetah is already connected';
else
    serverName = '192.168.3.100';
    disp(sprintf('Connecting to %s...', serverName));
    success = NlxConnectToServer(serverName);
    if success ~= 1
        disp(sprintf('FAILED connect to %s.', serverName));
        return;
    else
        disp(sprintf('Connected to Cheetah at %s.', serverName));
    end
end 
%% Read WAV file from filesystem:
wavfilename = 'SingleShutterClick_short.wav';
[y, freq] = psychwavread(wavfilename);
wavedata = y';
nrchannels = size(wavedata,1); % Number of rows == number of channels.

if nrchannels < 2
    wavedata = [wavedata ; wavedata];
    nrchannels = 2;
end

%%
InitializePsychSound;
device = [];
pahandle = PsychPortAudio('Open', device, [], 0, freq, nrchannels);
PsychPortAudio('FillBuffer', pahandle, wavedata);

%%
repetitions = 1;
t1 = PsychPortAudio('Start', pahandle, repetitions, 0, 1);

global t 
t = timer;
t.StartDelay = 3;
t.TimerFcn = @(~,~)disp('3 seconds have elapsed');
t.Period = 3;

t = timer('TimerFcn', 'stat=false; playShutterSoundScript_PTB',... 
                 'StartDelay',4);

