% JJS. 2022-11-15.

for iRep = 1:5   % first rep is just a 3 second wait
    if iRep == 1
        disp(num2str(iRep))
        playShutterSoundScript_PTB
    else
        disp(num2str(iRep))
        start(t)
        
        stat=true;
        while(stat==true)
            disp('.')
            pause(1)
        end
    end
end





