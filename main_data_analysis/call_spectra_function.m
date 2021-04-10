%% call spectra function

clearvars;
close all
clc

n_perm = 10000;

want_all_trls = false;
addpath('/home/elio/toolboxes/MATLAB/circstat-matlab/')

% frequency limits to apply multiple corrections
foi = [0, 15]; % whole range, for compatibility with previous version

rng(0)

%% 
if want_all_trls

    load('../HR_allsubj_allTrls.mat')
    load('../HR_allsubj_L02_alltrls.mat')

else
    load('../HR_allsubj.mat')
    load('../HR_allsubj_L02.mat')
end
n_subj = size(HR_mat,3);

load0 = squeeze(HR_mat(1,:,:));
load2 = squeeze(HR_mat(2,:,:));
load4 = squeeze(HR_mat(3,:,:));

%% create required strct

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


spctr_load0 = cmpt_beh_spectra(load0, params);
spctr_load2 = cmpt_beh_spectra(load2, params);
spctr_load4 = cmpt_beh_spectra(load4, params);

spctr_load02 = cmpt_beh_spectra(HR_L02, params);

mask_freq_interest = (spctr_load0.freqs >= min(foi)) & (spctr_load0.freqs <= max(foi));

%% find maxima frequencies


colors = [0 204 204;
         255 51 153]/255;




%% see what's going on

fin_L0 = mean(spctr_load0.spctr_out, 2);
fin_L2 = mean(spctr_load2.spctr_out, 2);
fin_L4 = mean(spctr_load4.spctr_out, 2);
fin_L02 = mean(spctr_load02.spctr_out, 2);

% and compare it with phase locked signals..
PL_L0 = abs(sum(spctr_load0.cmplx_out, 2))/n_subj;
PL_L2 = abs(sum(spctr_load2.cmplx_out, 2))/n_subj;
PL_L4 = abs(sum(spctr_load4.cmplx_out, 2))/n_subj;
PL_L02 = abs(sum(spctr_load02.cmplx_out, 2))/n_subj;


%% start permutations

loadconds = 2;

[big_mat_perm, phaselocked_perms_bigmat] = deal(nan([size(spctr_load0.spctr_out),...
    n_perm, loadconds]));
    
HR_mat_perm = permute_from_raw_lowhighload(want_all_trls, n_perm);

emp_multcomp = nan(n_perm, loadconds);
    
for iPerm = 1:n_perm

    for iLoad = 1:loadconds

        shuffled_HR_mat = squeeze(HR_mat_perm(iLoad,:,:,iPerm));
        curr_spctr = cmpt_beh_spectra(shuffled_HR_mat, params);        
        big_mat_perm(:,:,iPerm,iLoad) = curr_spctr.spctr_out;
        phaselocked_perms_bigmat(:,:,iPerm,iLoad) = curr_spctr.cmplx_out;
       
        foo = 1;
        
        vect_avgspctr = abs(sum(curr_spctr.cmplx_out(mask_freq_interest, :), 2));
        emp_multcomp(iPerm, iLoad) = max(vect_avgspctr);
        
    end

    if mod(iPerm, 100) ==0
        fprintf('\n %d permutations', iPerm)
    end

end
   
perm_PL_L02 = squeeze(abs(sum(phaselocked_perms_bigmat(:,:,:,1), 2))/n_subj);
perm_PL_L4 = squeeze(abs(sum(phaselocked_perms_bigmat(:,:,:,2), 2))/n_subj);

PL_L02_95 = prctile(perm_PL_L02, 95, 2);
PL_L4_95 = prctile(perm_PL_L4, 95, 2);


[mat_multcomp4(:, 1), mat_multcomp4(:, 2)] = ecdf(emp_multcomp(:, 2)/n_subj);

[mat_multcomp02(:, 1), mat_multcomp02(:, 2)] = ecdf(emp_multcomp(:, 1)/n_subj);





%%
set(groot, 'defaultAxesFontSize',14)


figure; 

subplot(2, 3, 1:2); hold on
plot(spctr_load0.freqs, PL_L02, 'Color', colors(1, :), 'LineWidth', 3)
plot(spctr_load0.freqs, PL_L02_95, 'r--', 'LineWidth', 2)
legend('low load [0 & 2]', '95 percentile permutations')
title('phase locked sum')
xlim(minmax(spctr_load0.freqs))
ylabel('amplitude (a.u.)')

subplot(2, 3, 4:5); hold on
plot(spctr_load0.freqs, PL_L4, 'Color', colors(2, :), 'LineWidth', 3)
plot(spctr_load0.freqs, PL_L4_95, 'r--', 'LineWidth', 2)
legend('high load [4]', '95 percentile permutations')
xlim(minmax(spctr_load0.freqs))
xlabel('frequency [Hz]')
ylabel('amplitude (a.u.)')

ax = subplot(2, 3, 3); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)
[peakfreqs_L02, cntr_freq, pval, Zcirc_L02] = peakfreqs(spctr_load02, pol_ax, ...
    mask_freq_interest, colors(1, :));


title(sprintf('Low load [0 | 2]\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq, pval), 'FontSize', 14)

ax = subplot(2, 3, 6); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)
[peakfreqs_L4, cntr_freq, pval, Zcirc_L4] = peakfreqs(spctr_load4, pol_ax, ...
    mask_freq_interest, colors(2, :));
title(sprintf('High Load [4]\n F=%0.2f Hz, Rayleigh p=%0.3f', cntr_freq, pval), 'FontSize', 14)



%% pvalues

pvals_L02 = compute_pvals_and_mc(PL_L02, ...
    squeeze(abs(sum(phaselocked_perms_bigmat(:, :, :, 1), 2))/n_subj), ...
    spctr_load0.freqs)

pvals_L4 = compute_pvals_and_mc(PL_L4, ...
    squeeze(abs(sum(phaselocked_perms_bigmat(:, :, :, 2), 2))/n_subj), ...
    spctr_load0.freqs)


%% cowans K and maxima
% load cowan's K 
load('../CowansK.mat')


[r_MAXcorr, p_MAXcorr] = corr(peakfreqs_L4, cowanK_mat(:, 2))
[r_LOWcorr, p_LOWcorr] = corr(peakfreqs_L02, cowanK_mat(:, 2))


figure; 
subplot(1, 2, 1)
scatter(peakfreqs_L4, cowanK_mat(:, 2))
subplot(1, 2, 2)
scatter(peakfreqs_L02, cowanK_mat(:, 2))




%% plot dynamic range after detrend 
% still keep the three load condition separate here

det_HR = nan(size(HR_mat));

str_load = {'load0', 'load2', 'load4'};

for iLoad = 1:3
    
    currMat = squeeze(HR_mat(iLoad,:,:));
    det_HR(iLoad,:,:) = do_detrend(currMat, size(HR_mat, 3), params);
    
end

min_det = squeeze(min(det_HR, [], 2));
max_det = squeeze(max(det_HR, [], 2));

figure();

mat_plotpos = [1, 2, 3; 5, 6, 7; 9, 10, 11];


for icond =1:3
    
    
    subplot(3, 4, mat_plotpos(icond, :))
    
    mat_bar = [min_det(icond, :); max_det(icond, :)]';
    bar(mat_bar, 'stacked')
    if icond==3
        xlabel('subjects')
    end
    ylabel('dynamic range')
    title(str_load{icond})
    ylim([-.4, .4])
    
    subplot(3, 4, icond*4)
    avg_bar = mean(mat_bar);
    err_bar = std(mat_bar) / sqrt(22);
    bar(avg_bar(1)); hold on
    bar(avg_bar(2))
    errorbar(avg_bar(1), err_bar(1), 'k', 'LineWidth', 2)
    errorbar(avg_bar(2), err_bar(2), 'k', 'LineWidth', 2)
    if icond == 1
        title('mean and std error')
    end
    xticks([])
    ylim([-.4, .4])
    
end









