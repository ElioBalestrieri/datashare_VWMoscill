%% call spectra function

clearvars;
close all
clc

n_perm = 1000;

permutation = 'raw';
want_all_trls = false;
addpath('/home/ebalestr/toolboxes/CircStat/')


%% 
if want_all_trls
    load('../HR_allsubj_allTrls.mat')
else
    load('../HR_allsubj.mat')
end
n_subj = size(HR_mat,3);

load0 = squeeze(HR_mat(1,:,:));
load2 = squeeze(HR_mat(2,:,:));
load4 = squeeze(HR_mat(3,:,:));

%% create required strct

params = [];
params.detrend_flag = 2;
params.window = 'hanning';
params.power = 0;
params.zero_pad = 7;
params.subj_dim = 2;
params.time_bins = (.15:.04:.75)';
params.f_sample = 25;
params.verbose = -1;
params.lp_filter = 0;


spctr_load0 = cmpt_beh_spectra(load0, params);
spctr_load2 = cmpt_beh_spectra(load2, params);
spctr_load4 = cmpt_beh_spectra(load4, params);



%% find maxima frequencies

figure()

ax = subplot(1, 3, 1); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)
[peakfreqs_L0, maxfreq, pval] = peakfreqs(spctr_load0, pol_ax);
title(sprintf('Load 0, F=%0.2f, Rayleigh p=%0.3f', maxfreq, pval))

ax = subplot(1, 3, 2); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)
[peakfreqs_L2, maxfreq, pval] = peakfreqs(spctr_load2, pol_ax);
title(sprintf('Load 2, F=%0.2f, Rayleigh p=%0.3f', maxfreq, pval))

ax = subplot(1, 3, 3); 
pol_ax = polaraxes('Units', ax.Units, 'Position',ax.Position);
delete(ax)
[peakfreqs_L4, maxfreq, pval] = peakfreqs(spctr_load4, pol_ax);
title(sprintf('Load 4, F=%0.2f, Rayleigh p=%0.3f', maxfreq, pval))

%%

condnames = {'load0', 'load2', 'load4'};

% load cowan's K 
load('../CowansK.mat')

% collect all peak frequencies in one mat for ANOVA
mat_peak_freqs = [peakfreqs_L0, peakfreqs_L2, peakfreqs_L4];

% correlate individual peaks with cowans K in respoective WM conditions
[rho2, corr_p2] = corr(peakfreqs_L2, cowanK_mat(:, 1));
[rho4, corr_p4] = corr(peakfreqs_L4, cowanK_mat(:, 2));

figure();
subplot(1, 2, 1)
scatter(peakfreqs_L2, cowanK_mat(:, 1))
ylabel('Cowan K')
xlabel('peak frequency')
title('load 2')
subplot(1, 2, 2)
scatter(peakfreqs_L4, cowanK_mat(:, 2))
xlabel('peak frequency')
title('load 4')

figure()
[outTableANOVA, ~, ~] = rm1W_ANOVA_adapted(mat_peak_freqs,...
    condnames,0,1,'peak frequencies')


% compute phase consistency for maxima




%% see what's going on

fin_L0 = mean(spctr_load0.spctr_out, 2);
fin_L2 = mean(spctr_load2.spctr_out, 2);
fin_L4 = mean(spctr_load4.spctr_out, 2);

% and compare it with phase locked signals..
PL_L0 = abs(sum(spctr_load0.cmplx_out, 2))/n_subj;
PL_L2 = abs(sum(spctr_load2.cmplx_out, 2))/n_subj;
PL_L4 = abs(sum(spctr_load4.cmplx_out, 2))/n_subj;

% find maxima 



figure()
subplot(2, 1, 1); hold on
plot(spctr_load0.freqs, fin_L0)
plot(spctr_load0.freqs, fin_L2)
plot(spctr_load0.freqs, fin_L4)

subplot(2, 1, 2); hold on
plot(spctr_load0.freqs, PL_L0)
plot(spctr_load0.freqs, PL_L2)
plot(spctr_load0.freqs, PL_L4)


[outANOVAS, structT] =  SPECTRAL_analysis(spctr_load0.spctr_out, spctr_load2.spctr_out,...
    spctr_load4.spctr_out, spctr_load0.freqs);

%% start permutations

[big_mat_perm, phaselocked_perms_bigmat] = deal(nan([size(spctr_load0.spctr_out),...
    n_perm, 3]));
    
HR_mat_perm = permute_from_raw(want_all_trls, n_perm);
    
for iPerm = 1:n_perm

    for iLoad = 1:3

        shuffled_HR_mat = squeeze(HR_mat_perm(iLoad,:,:,iPerm));
        curr_spctr = cmpt_beh_spectra(shuffled_HR_mat, params);        
        big_mat_perm(:,:,iPerm,iLoad) = curr_spctr.spctr_out;
        phaselocked_perms_bigmat(:,:,iPerm,iLoad) = curr_spctr.cmplx_out;
       
    end

    if mod(iPerm, 100) ==0
        fprintf('\n %d permutations', iPerm)
    end

end
   
perm_L0 = squeeze(mean(big_mat_perm(:,:,:,1), 2));
perm_L2 = squeeze(mean(big_mat_perm(:,:,:,2), 2));
perm_L4 = squeeze(mean(big_mat_perm(:,:,:,3), 2));

perm_PL_L0 = squeeze(abs(sum(phaselocked_perms_bigmat(:,:,:,1), 2))/n_subj);
perm_PL_L2 = squeeze(abs(sum(phaselocked_perms_bigmat(:,:,:,2), 2))/n_subj);
perm_PL_L4 = squeeze(abs(sum(phaselocked_perms_bigmat(:,:,:,3), 2))/n_subj);

L0_95 = prctile(perm_L0, 95, 2);
L2_95 = prctile(perm_L2, 95, 2);
L4_95 = prctile(perm_L4, 95, 2);

PL_L0_95 = prctile(perm_PL_L0, 95, 2);
PL_L2_95 = prctile(perm_PL_L2, 95, 2);
PL_L4_95 = prctile(perm_PL_L4, 95, 2);


%% plot 95 percentiles with phase locked sum 
figure()    

subplot(1, 3, 1); hold on
plot(spctr_load0.freqs, PL_L0, 'LineWidth', 2)
plot(spctr_load0.freqs, PL_L0_95, 'LineWidth', 2)

subplot(1, 3, 2); hold on
plot(spctr_load0.freqs, PL_L2, 'LineWidth', 2)
plot(spctr_load0.freqs, PL_L2_95, 'LineWidth', 2)

subplot(1, 3, 3); hold on
plot(spctr_load0.freqs, PL_L4, 'LineWidth', 2)
plot(spctr_load0.freqs, PL_L4_95, 'LineWidth', 2)


%% redo figure applying omnibus correction in the interval 2-10 Hz
% +mseb

limFreqsLgcl = spctr_load0.freqs>2 & spctr_load0.freqs<10;

L0_95(limFreqsLgcl==0)=nan;
L2_95(limFreqsLgcl==0)=nan;
L4_95(limFreqsLgcl==0)=nan;

mask.lglc.L0 = L0_95==max(L0_95);
mask.lglc.L2 = L2_95==max(L2_95);
mask.lglc.L4 = L4_95==max(L4_95);

mask.max.L0 = max(L0_95);
mask.max.L2 = max(L2_95);
mask.max.L4 = max(L4_95);


str_load = {'load 0', 'load 2', 'load 4'};
cols = [0 204 204;
        127 0 255;
        255 51 153]/255;

for iLoad = 1:3
    
    currMat = squeeze(HR_mat(iLoad,:,:));
    det_HR(iLoad,:,:) = apply_detrend(currMat, params);
    
end

%% plot dynamic range after detrend
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
    
    
end






spctr_mat = cat(3, spctr_load0.spctr_out, spctr_load2.spctr_out,spctr_load4.spctr_out)

%% get p vals

dists.L0 = perm_L0(mask.lglc.L0,:);
dists.L2 = perm_L2(mask.lglc.L2,:);
dists.L4 = perm_L4(mask.lglc.L4,:);

peaks = max([fin_L0, fin_L2, fin_L4])

fn = fieldnames(dists);
for iLoad = 1:3
    
    [arrProb(:,1,iLoad), arrProb(:,2,iLoad)] = ecdf(dists.(fn{iLoad}));
    
    idxVal = find(arrProb(:,2,iLoad)>peaks(iLoad),1)-1;
    pVals(iLoad) = 1-arrProb(idxVal, 1, iLoad);
    
end
    
    

%% nice plot
figure
cPlot = 0;
lineProps = [];
vect_thresh = [mask.max.L0, mask.max.L2, mask.max.L4];

for iPlot = 1:3
    
    PlotStr(iPlot).RTavg = squeeze(mean(HR_mat(iPlot,:,:),3)); % is Hit rate, not Reaction times :)
    PlotStr(iPlot).RTstdERR = squeeze(std(HR_mat(iPlot,:,:),[],3))/sqrt(n_subj);
    
    PlotStr(iPlot).detRTavg = squeeze(mean(det_HR(iPlot,:,:),3));
    PlotStr(iPlot).detRTstdERR = squeeze(std(det_HR(iPlot,:,:),[],3))/sqrt(n_subj);
    
    PlotStr(iPlot).spctrAVG = squeeze(mean(spctr_mat(:,:,iPlot),2));
    PlotStr(iPlot).spctrStdERR = squeeze(std(spctr_mat(:,:,iPlot),[],2))/sqrt(n_subj);
    
    cPlot = cPlot+1;
    ax = subplot(6,2,cPlot)
    lineProps.col = {cols(iPlot,:)};
    lineProps.style = '-';
    mseb(params.time_bins', PlotStr(iPlot).RTavg, PlotStr(iPlot).RTstdERR, lineProps, 1)
    xlim(minmax(params.time_bins'))
    title(str_load{iPlot})
    ylabel('HR')
    set(ax(1),'XTickLabel','')
    set(ax(1),'XColor',get(gca,'Color'))
    set(ax(1),'box','off')
    
    
    cPlot = cPlot+2;
    ax = subplot(6,2,cPlot); hold on;
    lineProps.style = '-';
    mseb(params.time_bins', PlotStr(iPlot).detRTavg, PlotStr(iPlot).detRTstdERR, lineProps, 1)
    plot(params.time_bins, zeros(numel(params.time_bins), 1), 'k', 'LineWidth',.5) 
    xlim(minmax(params.time_bins'))
    ylabel({'detrended', 'HR'})
    
    if iPlot<3
        set(ax(1),'XTickLabel','')
        set(ax(1),'XColor',get(gca,'Color'))
        
    end
    
    if iPlot == 3
        xlabel('time (s)')
    end
    
    ax = subplot(6,2,[cPlot-1,cPlot+1]);
    lineProps.style = '-';
    mseb(spctr_load0.freqs, PlotStr(iPlot).spctrAVG', PlotStr(iPlot).spctrStdERR', lineProps, 1)
    
    hold on    
    plot([0 12.5], [vect_thresh(iPlot),vect_thresh(iPlot)],'--k','LineWidth', 2) 
    
    xlim([2 10])
    ylim([.15 .35])
    
    if iPlot ==1
        title('spectra')
    end
    
    if iPlot<3
        set(ax(1),'XTickLabel','')
        set(ax(1),'XColor',get(gca,'Color'))
    end
        
    if iPlot == 3
        xlabel('frequency (Hz)')
    end
    ylabel('fft amplitude')

    cPlot = cPlot+1;
    
    
end
    

%% final figure -truth moment-

L0_95 = prctile(perm_L0, 95, 2);
L2_95 = prctile(perm_L2, 95, 2);
L4_95 = prctile(perm_L4, 95, 2);


figure
subplot(1,3,1)
plot(spctr_load0.freqs, fin_L0,'b', 'LineWidth', 2); hold on;
plot(spctr_load0.freqs, L0_95,'r', 'LineWidth', 2); hold on;
title('load0')
xlabel('freqs (Hz)')
xlim(minmax(spctr_load0.freqs));
ylabel('power')

subplot(1,3,2)
plot(spctr_load0.freqs, fin_L2,'b', 'LineWidth', 2); hold on;
plot(spctr_load0.freqs, L2_95,'r', 'LineWidth', 2); hold on;
title('load2')
xlabel('freqs (Hz)')
xlim(minmax(spctr_load0.freqs));


subplot(1,3,3)
plot(spctr_load0.freqs, fin_L4,'b', 'LineWidth', 2); hold on;
plot(spctr_load0.freqs, L4_95,'r', 'LineWidth', 2); hold on;
title('load4')
xlabel('freqs (Hz)')
xlim(minmax(spctr_load0.freqs));

%% check uncorrected pVal for 

bigPerm = cat(3, perm_L0, perm_L2, perm_L4);
arrProb_ALL = nan(n_perm+1, 2, numel(spctr_load0.freqs), 3);
mat_spctr = [fin_L0, fin_L2, fin_L4];
for iLoad = 1:3
    
    for iFreq = 1:numel(L0_95)
        
        [arrProb_ALL(:,1,iFreq,iLoad), arrProb_ALL(:,2,iFreq,iLoad)] =...
            ecdf(bigPerm(iFreq, :, iLoad));
    
        idxVal = find(arrProb_ALL(:,2,iFreq,iLoad)>mat_spctr(iFreq, iLoad),1)-1;
        pVals_all(iFreq,iLoad) = 1-arrProb_ALL(idxVal, 1,iFreq, iLoad);
        
    end
    
end
        
        
pVals_all = [spctr_load0.freqs', pVals_all];

%% 

