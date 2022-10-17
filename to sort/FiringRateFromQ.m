function [F] = FiringRateFromQ(Q, S, dt, varargin)
%2020-03-09. JJS. Create a ctsd of Firing rate from Q matrix. 
%   Detailed explanation goes here

windowSize = 1; % seconds
gaussBlur = 0.1; % seconds
normalize = 0; % 1 = yes, 0 = no;

extract_varargin;

nC = length(S);

% Q = MakeQfromS(S, dt, 'Tstart', tStart, 'Tend', tEnd);
% QD = Data(Q);
% QD = full(QD);

QD = Q.data;

x = -windowSize:dt:windowSize;
y = normpdf(x,0,gaussBlur)';
y = y/sum(y);

for iC = 1:nC
    QD(:,iC) = conv2(QD(:,iC), y, 'same'); 
end

if normalize == 1;
    for iC = 1:nC
        QD(:,iC) = QD(:,iC)/sum(QD(:,iC));
    end
end
        
F = tsd(Q.tvec, QD);



end

