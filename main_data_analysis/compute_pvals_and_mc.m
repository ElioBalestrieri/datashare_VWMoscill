function out_struct = compute_pvals_and_mc(data, perms, freqs, varargin)

if ~isempty(varargin) % give the possibility of select a specific index
    
    max_peak = data(varargin{1});

else
    
    % find absolute maximum
    max_peak = max(data);

end


mask_max = data == max_peak;


% select the permutations for the peak frequency
perms_freqspec = perms(mask_max, :)';

% generate empirical CDF
[mat_ecdf(:, 1), mat_ecdf(:, 2)] = ecdf(perms_freqspec);

% generate uncorrected pvalue
idx_sign = find(mat_ecdf(:, 2) <= max_peak, 1, 'last');
pval_raw = 1 - mat_ecdf(idx_sign, 1);

% foi 
mask_foi = freqs>2 & freqs<10;
num_comp_in_range = sum(mask_foi);

% local maxima
locmax = islocalmax([data; 0]);


out_struct.freq = freqs(mask_max);
out_struct.p_uncorrected = pval_raw;
out_struct.p_corr_bonf_range = pval_raw*num_comp_in_range;
out_struct.p_corr_bonf_peaknumbers = pval_raw*sum(locmax);
out_struct.p_corr_bonferroni = pval_raw*length(freqs);
out_struct.npeaks = sum(locmax);








end