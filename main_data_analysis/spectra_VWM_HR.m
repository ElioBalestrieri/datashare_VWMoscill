clear all

addpath('/home/elio/toolboxes/MATLAB/npy-matlab/npy-matlab')  

load('../HR_allsubj.mat')
VWM_trends_SOA = readNPY('../VWM_trends_SOA.npy');


%% 

nsubjs = size(HR_mat, 3);

params = [];
params.detrend_flag = 2;
params.window = [];
params.power = 0;
params.zero_pad = 2;
params.subj_dim = 2;
params.time_bins = (.15:.04:.75)';
params.f_sample = 25;
params.verbose = -1;
params.lp_filter = 0;


spctr_load0 = cmpt_beh_spectra(squeeze(HR_mat(1, :, : )), params);
spctr_load2 = cmpt_beh_spectra(squeeze(HR_mat(2, :, : )), params);
spctr_load4 = cmpt_beh_spectra(squeeze(HR_mat(3, :, : )), params);

spctr_VWM_load2 = cmpt_beh_spectra(squeeze(VWM_trends_SOA(:, :, 1)), params);
spctr_VWM_load4 = cmpt_beh_spectra(squeeze(VWM_trends_SOA(:, :, 2)), params);


cols = [0 204 204;
        127 0 255;
        255 51 153]/255;


avg_VWM_series{1} = abs(sum(spctr_VWM_load2.cmplx_out, 2))/nsubjs;

avg_VWM_series{2} = abs(sum(spctr_VWM_load4.cmplx_out, 2))/nsubjs;




%% permute

n_iter = 10000;

VWM_permuted = VWM_permute_from_raw(n_iter);

% compute permuted spectra
[permspectra_L2,permspectra_L4] = deal(nan(length(avg_VWM_series{1}), ...
                                       n_iter));

for iperm = 1:n_iter
    
    L2_spctr = cmpt_beh_spectra(squeeze(VWM_permuted(1, :, :, iperm)), params);
    L4_spctr = cmpt_beh_spectra(squeeze(VWM_permuted(2, :, :, iperm)), params);
    
    permspectra_L2(:, iperm) = abs(sum(L2_spctr.cmplx_out, 2))/nsubjs;
    permspectra_L4(:, iperm) = abs(sum(L4_spctr.cmplx_out, 2))/nsubjs;    
        
end

thresh_load2 = prctile(permspectra_L2, 95, 2);
thresh_load4 = prctile(permspectra_L4, 95, 2);

%%

figure; 

subplot(1, 2, 1); hold on
plot(spctr_load0.freqs, avg_VWM_series{1}, 'LineWidth', 3)
plot(spctr_load0.freqs, thresh_load2, 'LineWidth', 3)
title('VWM accuracy spectrum, load2')
xlabel('freq (Hz)')
ylabel('amplitude (a.u)')
xlim(minmax(spctr_load0.freqs))
legend('signal', 'threshold (95° percentile permutations')

subplot(1, 2, 2); hold on
plot(spctr_load0.freqs, avg_VWM_series{2}, 'LineWidth', 3)
plot(spctr_load0.freqs, thresh_load4, 'LineWidth', 3)
title('VWM accuracy spectrum, load4')
xlabel('freq (Hz)')
ylabel('amplitude (a.u)')
xlim(minmax(spctr_load0.freqs))





