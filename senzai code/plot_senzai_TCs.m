% Mapstouse = []; 

% figure
hold on
for iCell = 1:length(Mapstouse)
    plot(10:10:360, Mapstouse{iCell,1}.rate)
%     pause
%     disp('press any key to continue')
end

sig77 = pval_ppln_77 < .05; num_sig77 = sum(sig77); fract_sig_77 = num_sig77/length(Maps77); 
sig79 = pval_ppln_79 < .05; num_sig79 = sum(sig79); fract_sig_79 = num_sig79/length(Maps79); 
sig83 = pval_ppln_83 < .05; num_sig83 = sum(sig83); fract_sig_83 = num_sig83/length(Maps83); 
sig85 = pval_ppln_85 < .05; num_sig85 = sum(sig85); fract_sig_85 = num_sig85/length(Maps85); 
sig96 = pval_ppln_96 < .05; num_sig96 = sum(sig96); fract_sig_96 = num_sig96/length(Maps96); 
sig98 = pval_ppln_98 < .05; num_sig98 = sum(sig98); fract_sig_98 = num_sig98/length(Maps98); 

open_field_sig_fraction = [fract_sig_77 fract_sig_79 fract_sig_83 fract_sig_85 fract_sig_96 fract_sig_98];