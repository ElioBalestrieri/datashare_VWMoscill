function peakfreqs_vect = peakfreqs(spctr_loadx)

nsubjs = size(spctr_loadx.spctr_out, 2);
max_eachsubj = max(spctr_loadx.spctr_out);
peakfreqs_vect = nan(nsubjs, 1);

for isubj = 1:nsubjs
    
    this_mag = spctr_loadx.spctr_out(:, isubj);
    mask_ = this_mag==max_eachsubj(isubj);
    peakfreqs_vect(isubj) = spctr_loadx.freqs(mask_);

end





end