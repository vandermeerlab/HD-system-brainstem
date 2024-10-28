figure
hold on
for iCell = 1:length(Maps)
    plot(Maps{iCell,1}.rate)
    disp('press any key to continue')
end


for iCell = 1:length(Mapstouse)
    plot(10:10:360, Mapstouse{iCell,1}.rate)
%     pause
%     disp('press any key to continue')
end