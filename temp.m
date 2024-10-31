plot(M77_AHVtsd.tvec, M77_AHVtsd.data, '.'); hold on
plot(M77_heading_degreesTSD_lowAHV.tvec, M77_heading_degreesTSD_lowAHV.data, 'r.')
plot(M77_saccade.t, M77_heading_degreesTSD_lowAHV.data(M77_saccade.t))