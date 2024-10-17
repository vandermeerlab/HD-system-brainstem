figure
hold on
for iCell = 1:length(Maps)
    plot(Maps{iCell,1}.rate)
    disp('press any key to continue')
end