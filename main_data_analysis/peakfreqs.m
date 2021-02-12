function [peakfreqs_vect, maxfreq, pval_circ] = peakfreqs(spctr_loadx, pol_ax)

nsubjs = size(spctr_loadx.spctr_out, 2);
max_eachsubj = max(spctr_loadx.spctr_out);
peakfreqs_vect = nan(nsubjs, 1);

for isubj = 1:nsubjs
    
    this_mag = spctr_loadx.spctr_out(:, isubj);
    mask_ = this_mag==max_eachsubj(isubj);
    peakfreqs_vect(isubj) = spctr_loadx.freqs(mask_);

end


% compute phase consistency at population maxima
avg_ampl = mean(spctr_loadx.spctr_out, 2);
% avg_ampl = abs(sum(spctr_loadx.cmplx_out, 2));


max_mask = avg_ampl==max(avg_ampl);
maxfreq = spctr_loadx.freqs(max_mask);

angles_ = angle(spctr_loadx.cmplx_out(max_mask, :));
pval_circ = circ_rtest(angles_);


out_pol = polarhistogram(pol_ax, angles_, 10,...
    'FaceColor',  [.1 .1 .1],...
    'Normalization', 'probability');

hold on;
angle_avg = circ_mean(angles_');

polarscatter(angle_avg, max(out_pol.Values), 30, [0, 153, 153]/255, 'filled')
polarplot([0, angle_avg], [0,  max(out_pol.Values)], 'LineWidth', 2, ...
    'Color', [0, 153, 153]/255)





end