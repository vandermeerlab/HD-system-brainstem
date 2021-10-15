
box off
set(gca,'FontSize',16,'YTick',[],'XTick',[])
ylabel('Trial number')
title('Saccade PETH')
hold on
for iLev = 1:2    
        s_idx = find(S.t{iC} > task.leverEvt.t{iLev}(iT) -2 & S.t{iC} < task.leverEvt.t{iLev}(iT) + 15);        
        if ~isempty(s_idx)            
            spikes = S.t{iC}(s_idx) - task.leverEvt.t{iLev}(iT);            
            if length(spikes) == 2                
                plot([spikes(1) spikes(1)],[line_h line_h+2.8],'color',colors{iLev})
                plot([spikes(2) spikes(2)],[line_h line_h+2.8],'color',colors{iLev})            
            else
                plot([spikes spikes],[line_h line_h+2.8],'color',colors{iLev})            
            end
        end
        line_h = line_h + 3;    
    end
    line_h = line_h + 4;
end
xlim([-2 15]); ylim([0 line_h-3]);
plot([0 0],[0 line_h-3],'--r')