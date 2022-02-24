function out = encoder(a,b)

% function encoder - encodes the movements of quadrature encoder (HEDR-5420-ES214) and return
% a vector containing the direction of movements. 
% Input : 2 vectors - a and b, which are the direct output of the quadrature
% encoder and have the same size.
% Output : Ternary vector same size as input vectors. 0 - no movement. 1-forward. -1 backward.


out = zeros(1,length(a)); % the output vector
aup =  (a > 20); % a at up
achange = findChangeNoDC(aup,0.5); % points of state change
bup =  (b > 20); % b at up
bchange = findChangeNoDC(bup,0.5); % points of state change
bothup = and(aup,bup); % logical - both signals at up state
statechange = findChangeNoDC(bothup,0.5);

if length(statechange) > 1 % proceed only if there is movement (state change in encoder)
    % testing for direction change - go over all points and check which
    % channel preceeds the other.
    lastind = 0;
    count = 0;
    for states = 1:length(statechange) % at statechange both channels are at UP state
        % Negative direction
        if sum(aup((statechange(states)-1):((statechange(states)+1)))) < 3 % if channel 'a' was changed look for change in 'b' at previous time
            % look for previous positive change in 'b'
            bphase = find(and(bchange < statechange(states), bchange > lastind));
            if bphase > 0
                count = count - 1;
                out(statechange(states)) = -1;
            end
            
        end
        
        %Positive direction - A preceedes B
        if sum(bup((statechange(states)-1):((statechange(states)+1)))) < 3 % if channel 'b' was changed look for change in 'a' at previous time
            % look for previous positive change in 'a'
            aphase = find(and(achange < statechange(states), achange > lastind));
            if length(aphase) > 0
                count = count + 1;
                out(statechange(states)) = +1;
            end
        end

        lastind = statechange(states); % the last indice to search up to
    end
end % go into loop only if there is change in signal

% plots the output over the input vectors
figure(43); clf;
plot(out*20,'r*') ; hold on; plot(a) ; plot(b)
%

end



