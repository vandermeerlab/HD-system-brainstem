    function simplePETH
    % bar graph
    m = histc(outputS, outputIT);
    bar(outputIT,m/cfg.dt/length(t));
    % 	x = outputIT;
    % 	m = nanmean(1./outputID);
    % 	se =  nanstd(1./outputID)/sqrt(nT+1);
    % 	plot(x,m,'b',x,m+se,'r:',x,m-se,'r:');
    set(gca, 'XLim', cfg.window);
    ylabel('FR (Hz)')
    xlabel('peri-event (sec)');% mean frequency line
    end