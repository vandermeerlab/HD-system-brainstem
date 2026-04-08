% Switch from unix start time to session start time
events = LoadEvents([]);
Events = [];
for iE = 1:length(events.t)
    Events.t{iE} = events.t{iE} - events.t{5};  % start record should be event.t{5}
end
Events.type = events.type;
Events.label = events.label; 


MarkerSize = 30; 
clf
plot(Events.t{2}, 1:length(Events.t{2}), '.', 'Color', 'g', 'MarkerSize', MarkerSize)  % Airpuff tone start
hold on
plot(Events.t{8}, 1:length(Events.t{8}), '.', 'Color', 'k', 'MarkerSize', MarkerSize, 'Marker', '|') % Airpuff tone end   TTL Output on AcqSystem1_0 board 0 port 0 value (0x0000).

plot(Events.t{1}, 1:length(Events.t{1}), '.', 'Color', 'r', 'MarkerSize', MarkerSize)  % Airpuff delivery 

plot(Events.t{4}, 1:length(Events.t{4}), '.', 'Color', 'c', 'MarkerSize', MarkerSize) % Control Tone Start

plot(Events.t{9}, 1:length(Events.t{9}), '.', 'Color', 'k', 'MarkerSize', MarkerSize, 'Marker', '|') % Control Tone End

plot(Events.t{3}, 1:length(Events.t{3}), '.', 'Color', 'm', 'MarkerSize', MarkerSize) % Control puff delivered

legend('Airpuff tone start', 'Airpuff tone end', 'Airpuff delivery', 'Control Tone start', 'Control Tone End', 'Control puff delivery')

set(gca, 'FontSize', 22)












% for iSess = 1:length(fd)
%     [~, b, ~] = fileparts(fd{iSess});
%     fdr{iSess} = b;
% end

for iCell = 1:length(Z)
    percent_n(iCell) = Z{iCell}.percent_n(1,end);
    percent_t(iCell) = Z{iCell}.percent_t(1,end);
end