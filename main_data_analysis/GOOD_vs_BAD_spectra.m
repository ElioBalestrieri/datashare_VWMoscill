%% compute spectra for bad and good

clearvars
close all
clc

want_all_trls = false;

n_iter = 10000;

addpath('/home/elio/toolboxes/MATLAB/circstat-matlab/')

rng(0)

load('../CowansK.mat')

if want_all_trls

    load('../HR_allsubj_allTrls.mat')
    load('../HR_allsubj_L02_alltrls.mat')
    
else
    
    load('../HR_allsubj.mat')
    load('../HR_allsubj_L02.mat')
    
end
    
    
medianK = median(cowanK_mat(:,2));

lgcl_good = cowanK_mat(:,2)>=medianK;
lgcl_poor = cowanK_mat(:,2)<medianK;

HR_good_L0 = squeeze(HR_mat(1,:,lgcl_good));
HR_good_L2 = squeeze(HR_mat(2,:,lgcl_good));
HR_good_L4 = squeeze(HR_mat(3,:,lgcl_good));

HR_poor_L0 = squeeze(HR_mat(1,:,lgcl_poor));
HR_poor_L2 = squeeze(HR_mat(2,:,lgcl_poor));
HR_poor_L4 = squeeze(HR_mat(3,:,lgcl_poor));

HR_poor_L02 = HR_L02(:, lgcl_poor);
HR_good_L02 = HR_L02(:, lgcl_good);

%% compute spectra

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

% good spectra
spctr_good_L0 = cmpt_beh_spectra(HR_good_L0, params);
spctr_good_L2 = cmpt_beh_spectra(HR_good_L2, params);
spctr_good_L4 = cmpt_beh_spectra(HR_good_L4, params);
spctr_good_L02 = cmpt_beh_spectra(HR_good_L02, params);

% poor spectra
spctr_poor_L0 = cmpt_beh_spectra(HR_poor_L0, params);
spctr_poor_L2 = cmpt_beh_spectra(HR_poor_L2, params);
spctr_poor_L4 = cmpt_beh_spectra(HR_poor_L4, params);
spctr_poor_L02 = cmpt_beh_spectra(HR_poor_L02, params);

% first visualization
nsubj = sum(lgcl_poor); % same numerosity
avg_poor4 = abs(sum(spctr_poor_L4.cmplx_out, 2))/nsubj;
avg_good4 = abs(sum(spctr_good_L4.cmplx_out, 2))/nsubj;

avg_poor02 = abs(sum(spctr_poor_L02.cmplx_out, 2))/nsubj;
avg_good02 = abs(sum(spctr_good_L02.cmplx_out, 2))/nsubj;

foi = [0, 9];
mask_freq_interest = (spctr_poor_L0.freqs >= min(foi)) & (spctr_poor_L0.freqs <= max(foi));

%%



colors = [1, 0, 0;
          0, 1, 0];


figure;

set(groot, 'defaultAxesFontSize',14)

ax = subplot(2, 4, 1); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)

[peakfreqs_poorL02, cntr_freq, pval, Z] = peakfreqs(spctr_poor_L02, pol_ax, ...
    mask_freq_interest, colors(1, :));
title(sprintf('poor Low load [0 & 2]\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq, pval), 'FontSize', 14)

ax = subplot(2, 4, 4); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)

[peakfreqs_goodL02, cntr_freq, pval, Z] = peakfreqs(spctr_good_L02, pol_ax, ...
    mask_freq_interest, colors(2, :));
title(sprintf('good Low load [0 & 2]\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq, pval), 'FontSize', 14)


ax = subplot(2, 4, 5); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)

[peakfreqs_poorL02, cntr_freq, pval, Z] = peakfreqs(spctr_poor_L4, pol_ax, ...
    mask_freq_interest, colors(1, :));
title(sprintf('poor high load [4]\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq, pval), 'FontSize', 14)

ax = subplot(2, 4, 8); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)

[peakfreqs_goodL02, cntr_freq, pval, Z] = peakfreqs(spctr_good_L4, pol_ax, ...
    mask_freq_interest, colors(2, :));
title(sprintf('good high load [4]\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq, pval), 'FontSize', 14)



%% permutations
HR_mat_perm = permute_from_raw_lowhighload(want_all_trls, n_iter);

[perms_G, perms_P] = deal(nan(n_iter, length(spctr_good_L02.freqs), 2));

[mult_comp_G, mult_comp_P] = deal(nan(n_iter, 2));


for iperm = 1:n_iter
    
    for iload = 1:2
        
        this_good = squeeze(HR_mat_perm(iload, :, lgcl_good, iperm));      
        permspctr_good = cmpt_beh_spectra(this_good, params);
        temp_G = abs(sum(permspctr_good.cmplx_out, 2))/nsubj;
        perms_G(iperm, :, iload) = temp_G;
        
        mult_comp_G(iperm, iload) = max(temp_G(mask_freq_interest));

        this_poor = squeeze(HR_mat_perm(iload, :, lgcl_poor, iperm));      
        permspctr_poor = cmpt_beh_spectra(this_poor, params);
        temp_P = abs(sum(permspctr_poor.cmplx_out, 2))/nsubj;
        perms_P(iperm, :, iload) = temp_P;
        
        mult_comp_P(iperm, iload) = max(temp_P(mask_freq_interest));
        
    end
    
end

%% plots

upper_good_L02 = prctile(squeeze(perms_G(:, :, 1)), 95);
upper_poor_L02 = prctile(squeeze(perms_P(:, :, 1)), 95);
upper_good_L4 = prctile(squeeze(perms_G(:, :, 2)), 95);
upper_poor_L4 = prctile(squeeze(perms_P(:, :, 2)), 95);



subplot(2, 4, 2:3); hold on
plot(spctr_good_L0.freqs, avg_poor02, 'r', 'LineWidth', 3)
plot(spctr_good_L0.freqs, upper_poor_L02, '--r', 'LineWidth', 3)

plot(spctr_good_L0.freqs, avg_good02, 'g', 'LineWidth', 3)
plot(spctr_good_L0.freqs, upper_good_L02, '--g', 'LineWidth', 3)

legend('poor performers (PP)', 'PP sign. thresh.', ...
       'good performers (GP)', 'GP sign. thresh.')

title('Low load [0 or 2]')
ylabel('amplitude (a.u.)')
xlim(minmax(spctr_good_L0.freqs))

subplot(2, 4, 6:7); hold on
plot(spctr_good_L0.freqs, avg_poor4, 'r', 'LineWidth', 3)
plot(spctr_good_L0.freqs, upper_poor_L4, '--r', 'LineWidth', 3)

plot(spctr_good_L0.freqs, avg_good4, 'g', 'LineWidth', 3)
plot(spctr_good_L0.freqs, upper_good_L4, '--g', 'LineWidth', 3)

title('High Load [4]')

xlabel('freqs (Hz)')
ylabel('amplitude (a.u.)')
xlim(minmax(spctr_good_L0.freqs))


%% multiple comparisons

pvals_L02_poor = compute_pvals_and_mc(avg_poor02, ...
    squeeze(perms_P(:, :, 1))', ...
    spctr_good_L0.freqs, 7)

pvals_L02_poor_nyq = compute_pvals_and_mc(avg_poor02, ...
    squeeze(perms_P(:, :, 1))', ...
    spctr_good_L0.freqs)

pvals_L02_good = compute_pvals_and_mc(avg_good02, ...
    squeeze(perms_G(:, :, 1))', ...
    spctr_good_L0.freqs)

pvals_L4_poor = compute_pvals_and_mc(avg_poor4, ...
    squeeze(perms_P(:, :, 2))', ...
    spctr_good_L0.freqs)

pvals_L4_good = compute_pvals_and_mc(avg_good4, ...
    squeeze(perms_G(:, :, 2))', ...
    spctr_good_L0.freqs)


