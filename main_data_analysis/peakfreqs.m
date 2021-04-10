function [peakfreqs_vect, which_freq, pval_circ, Z_circ] = peakfreqs(spctr_loadx, ...
    pol_ax, mask_freq_interest, facecolor)

nsubjs = size(spctr_loadx.spctr_out, 2);
max_eachsubj = max(spctr_loadx.spctr_out(mask_freq_interest, :));
peakfreqs_vect = nan(nsubjs, 1);
redfreqs = spctr_loadx.freqs(mask_freq_interest);

for isubj = 1:nsubjs
    
    this_mag = spctr_loadx.spctr_out(mask_freq_interest, isubj);
    mask_ = this_mag==max_eachsubj(isubj);
    peakfreqs_vect(isubj) = redfreqs(mask_);

end


% compute phase consistency at population maxima
avg_ampl = abs(sum(spctr_loadx.cmplx_out, 2))/nsubjs;
max_mask = avg_ampl==max(avg_ampl(mask_freq_interest));

which_freq = spctr_loadx.freqs(max_mask);

angles_ = angle(spctr_loadx.cmplx_out(max_mask, :));
[pval_circ, Z_circ] = circ_rtest(angles_);

out_pol = polarhistogram(pol_ax, angles_, 10,...
    'FaceColor',  facecolor,...
    'Normalization', 'probability');

hold on;
angle_avg = circ_mean(angles_');

polarscatter(angle_avg, max(out_pol.Values), 30, [0, 0, 0]/255, 'filled')
polarplot([0, angle_avg], [0,  max(out_pol.Values)], 'LineWidth', 2, ...
    'Color', [0, 0, 0]/255)


%% select phase in specific phase range across participants

% foo = 1;
% 
% freq_mask = abs(spctr_loadx.freqs-maxfreq)<2;
% freqs_selected = spctr_loadx.freqs(freq_mask);
% 
% 
% store_localpeaks = nan(nsubjs, 2);
% for isubj = 1:nsubjs
%     
%     % find max in range
%     vect_amp = spctr_loadx.spctr_out(freq_mask, isubj);
%     vect_cmplx = spctr_loadx.cmplx_out(freq_mask, isubj);
%     
%     max_cmplx = vect_cmplx(vect_amp == max(vect_amp));
%     
%     store_localpeaks(isubj, 1) = freqs_selected(vect_amp == max(vect_amp));
%     store_localpeaks(isubj, 2) = angle(max_cmplx);
%     
%     
%     
% end
% 
% 
% 
% %%
% 
% pval_circ = circ_rtest(store_localpeaks(:, 2));
% which_freq = mean(store_localpeaks(:, 1));
% 
% out_pol = polarhistogram(pol_ax, store_localpeaks(:, 2), 10,...
%     'FaceColor',  [.1 .1 .1],...
%     'Normalization', 'probability');
% 
% hold on;
% angle_avg = circ_mean(store_localpeaks(:, 2));
% 
% polarscatter(angle_avg, max(out_pol.Values), 30, [0, 153, 153]/255, 'filled')
% polarplot([0, angle_avg], [0,  max(out_pol.Values)], 'LineWidth', 2, ...
%     'Color', [0, 153, 153]/255)



end