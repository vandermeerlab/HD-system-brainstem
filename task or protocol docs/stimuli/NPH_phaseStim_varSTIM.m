

%% Script to send 100 stims pre-recording with 1 sec ISI
% Setup channel 1
chan_id = 1;
nstim = 100;
ConfirmBit = ProgramPulsePalParam(chan_id, 1, 0);  % parameter 1 is monophasic
ConfirmBit = ProgramPulsePalParam(chan_id, 4, 0.001);  % parameter 4 is pulse width
% Setup train with variable ISI between 0.5-3.5
V = 5*ones(1,nstim); % V and T need to be the same length. 
T = 3:1:102;
ConfirmBit = SendCustomPulseTrain(1,T,V); % first value can be 1 or 2. Can save two pulse trains. 
ConfirmBit = ProgramPulsePalParam(chan_id, 14, 1);  % parameter 14 is which custom channel to use 


%% Script to send 100 stims with variable ISI
% Setup channel 2
chan_id = 2;
nstim = 10;
ConfirmBit = ProgramPulsePalParam(chan_id, 1, 0);
ConfirmBit = ProgramPulsePalParam(chan_id, 4, 0.001);
% Setup train with variable ISI between 0.5-3.5
V = 5*ones(1,nstim);
var_ISI = 0.5 + 3.*rand(1,nstim);
var_ISI = round(var_ISI, 3);
T = ones(1,nstim)*4;
for i = 2:nstim
    T(i) = T(i-1) + var_ISI(i-1);
end
ConfirmBit = SendCustomPulseTrain(1,T,V);
ConfirmBit = ProgramPulsePalParam(chan_id, 14, 1);

%% Script to send 100 stims post-recording with 1 sec ISI
% Setup channel 3
chan_id = 3;
nstim = 100;
ConfirmBit = ProgramPulsePalParam(chan_id, 1, 0);
ConfirmBit = ProgramPulsePalParam(chan_id, 4, 0.001);
% Setup train with variable ISI between 0.5-3.5
V = 5*ones(1,nstim);
T = 3:1:102;
ConfirmBit = SendCustomPulseTrain(1,T,V);
ConfirmBit = ProgramPulsePalParam(chan_id, 14, 1);

%% Script to send 25 long stims podt-recording with 1 sec ISI
% Setup channel 4
chan_id = 4;
nstim = 25;
ConfirmBit = ProgramPulsePalParam(chan_id, 1, 0);
ConfirmBit = ProgramPulsePalParam(chan_id, 4, 0.05);
% Setup train with variable ISI between 0.5-3.5
V = 5*ones(1,nstim);
T = 3:1:27;
ConfirmBit = SendCustomPulseTrain(1,T,V);
ConfirmBit = ProgramPulsePalParam(chan_id, 14, 1);
