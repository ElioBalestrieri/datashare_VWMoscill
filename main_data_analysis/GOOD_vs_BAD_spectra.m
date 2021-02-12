%% compute spectra for bad and good

clearvars
% close all
clc

want_all_trls = false;




load('../CowansK.mat')

if want_all_trls

    load('../HR_allsubj_allTrls.mat')

else
    
    load('../HR_allsubj.mat')
    
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


%% compute spectra

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

% good spectra
spctr_good_L0 = cmpt_beh_spectra(HR_good_L0, params);
spctr_good_L2 = cmpt_beh_spectra(HR_good_L2, params);
spctr_good_L4 = cmpt_beh_spectra(HR_good_L4, params);

% poor spectra
spctr_poor_L0 = cmpt_beh_spectra(HR_poor_L0, params);
spctr_poor_L2 = cmpt_beh_spectra(HR_poor_L2, params);
spctr_poor_L4 = cmpt_beh_spectra(HR_poor_L4, params);

% average spectra AND stderr
% good
AVG_spctr_good_L0 = mean(spctr_good_L0.spctr_out,2);
stdERR_spctr_good_L0 = std(spctr_good_L0.spctr_out,[],2)/sqrt(size(spctr_good_L0.spctr_out,2));


AVG_spctr_good_L2 = mean(spctr_good_L2.spctr_out,2);
stdERR_spctr_good_L2 = std(spctr_good_L2.spctr_out,[],2)/sqrt(size(spctr_good_L2.spctr_out,2));

AVG_spctr_good_L4 = mean(spctr_good_L4.spctr_out,2);
stdERR_spctr_good_L4 = std(spctr_good_L4.spctr_out,[],2)/sqrt(size(spctr_good_L4.spctr_out,2));

% poor
AVG_spctr_poor_L0 = mean(spctr_poor_L0.spctr_out,2);
stdERR_spctr_poor_L0 = std(spctr_poor_L0.spctr_out,[],2)/sqrt(size(spctr_poor_L0.spctr_out,2));

AVG_spctr_poor_L2 = mean(spctr_poor_L2.spctr_out,2);
stdERR_spctr_poor_L2 = std(spctr_poor_L2.spctr_out,[],2)/sqrt(size(spctr_poor_L2.spctr_out,2));

AVG_spctr_poor_L4 = mean(spctr_poor_L4.spctr_out,2);
stdERR_spctr_poor_L4 = std(spctr_poor_L4.spctr_out,[],2)/sqrt(size(spctr_poor_L4.spctr_out,2));

%% permutations
n_subj = size(HR_mat,3);
n_iter = 1000;
HR_mat_perm = permute_from_raw(want_all_trls, n_iter);

perm_spectra_good(n_iter,3).foo = [];
perm_spectra_poor(n_iter,3).foo = [];

[perm_AND_det.poor, perm_AND_det.good] = deal(nan(16, 11,3, n_iter));



for iPerm = 1:n_iter
   
    for iLoad = 1:3
        current_good = squeeze(HR_mat_perm(iLoad, :, lgcl_good, iPerm));
        current_poor = squeeze(HR_mat_perm(iLoad, :, lgcl_poor, iPerm));
        
        perm_AND_det.good(:,:,iLoad, iPerm) = apply_detrend(current_good,params);
        perm_AND_det.poor(:,:,iLoad, iPerm) = apply_detrend(current_poor,params);
        
    
        current_spectra_good = cmpt_beh_spectra(current_good, params);
        perm_spectra_good(iPerm, iLoad).spctr_out = current_spectra_good.spctr_out;
          
        current_spectra_poor = cmpt_beh_spectra(current_poor, params);
        perm_spectra_poor(iPerm, iLoad).spctr_out = current_spectra_poor.spctr_out;
        
    end
    
    if mod(iPerm,100)==0
        fprintf('\n %d permutations', iPerm)
    end       
    
end

big_mat_good_L0 = squeeze(mean(cat(3,perm_spectra_good(:,1).spctr_out),2));
big_mat_good_L2 = squeeze(mean(cat(3,perm_spectra_good(:,2).spctr_out),2));
big_mat_good_L4 = squeeze(mean(cat(3,perm_spectra_good(:,3).spctr_out),2));

big_mat_poor_L0 = squeeze(mean(cat(3,perm_spectra_poor(:,1).spctr_out),2));
big_mat_poor_L2 = squeeze(mean(cat(3,perm_spectra_poor(:,2).spctr_out),2));
big_mat_poor_L4 = squeeze(mean(cat(3,perm_spectra_poor(:,3).spctr_out),2));

upper_good_L0 = prctile(big_mat_good_L0, 95, 2);
upper_good_L2 = prctile(big_mat_good_L2, 95, 2);
upper_good_L4 = prctile(big_mat_good_L4, 95, 2);


upper_poor_L0 = prctile(big_mat_poor_L0, 95, 2);
upper_poor_L2 = prctile(big_mat_poor_L2, 95, 2);
upper_poor_L4 = prctile(big_mat_poor_L4, 95, 2);

%% compute pVals

inputP.good.L0 = big_mat_good_L0;
inputP.good.L2 = big_mat_good_L2;
inputP.good.L4 = big_mat_good_L4;


inputP.poor.L0 = big_mat_poor_L0;
inputP.poor.L2 = big_mat_poor_L2;
inputP.poor.L4 = big_mat_poor_L4;

% put real values into the same structure
inputRealSpectra.good.L0 = AVG_spctr_good_L0; 
inputRealSpectra.good.L2 = AVG_spctr_good_L2;
inputRealSpectra.good.L4 = AVG_spctr_good_L4;

inputRealSpectra.poor.L0 = AVG_spctr_poor_L0; 
inputRealSpectra.poor.L2 = AVG_spctr_poor_L2;
inputRealSpectra.poor.L4 = AVG_spctr_poor_L4;

% get the highest 95 percentile theshold over interval 2-10 Hz
omnibInput.good.L0 = upper_good_L0;
omnibInput.good.L2 = upper_good_L2;
omnibInput.good.L4 = upper_good_L4;

omnibInput.poor.L0 = upper_poor_L0;
omnibInput.poor.L2 = upper_poor_L2;
omnibInput.poor.L4 = upper_poor_L4;


fn1 = fieldnames(inputP);
fn2 = fieldnames(inputP.good);


pValsGoodVsPoor = nan(numel(spctr_good_L0.freqs), 12);
c = 0;
for iSample = 1:2
    
    for iLoad = 1:3
        
        c = c+1;
        for iFreq = 1:numel(AVG_spctr_good_L0)
            
            dataVect = inputP.(fn1{iSample}).(fn2{iLoad})(iFreq,:);
            currVal = inputRealSpectra.(fn1{iSample}).(fn2{iLoad})(iFreq);
            [cdfMat(:,1), cdfMat(:,2)] = ecdf(dataVect);
            idxVal = find(cdfMat(:,2)>currVal,1)-1;
            
            if idxVal==0
                idxVal = 1;
            elseif isempty(idxVal)
                idxVal = n_iter;
            end
            
            
            pValsGoodVsPoor(iFreq, c) = 1-cdfMat(idxVal,1);
            
            % do the same thing but with the highest 95 percentile
            % threshold
            
            if spctr_good_L0.freqs(iFreq)>2 && spctr_good_L0.freqs(iFreq)<10
                lgcl_freqs = spctr_good_L0.freqs>2 & spctr_good_L0.freqs<10;
                
                thisPrctiles = omnibInput.(fn1{iSample}).(fn2{iLoad});
                currMax = max(thisPrctiles(lgcl_freqs));
                max_perms = find(thisPrctiles==currMax);
                
                dataVect = inputP.(fn1{iSample}).(fn2{iLoad})(max_perms,:);
                currVal = inputRealSpectra.(fn1{iSample}).(fn2{iLoad})(iFreq);
                [cdfMat(:,1), cdfMat(:,2)] = ecdf(dataVect);
                idxVal = find(cdfMat(:,2)>currVal,1)-1;

                if idxVal==0
                    idxVal = 1;
                elseif isempty(idxVal)
                    idxVal = n_iter;
                end


                pValsGoodVsPoor(iFreq, c+6) = 1-cdfMat(idxVal,1);

            end
            
        end
        
    end
    
    
end


pValsGoodVsPoor = [spctr_good_L0.freqs', pValsGoodVsPoor];

pValsGoodVsPoor = array2table(pValsGoodVsPoor,'VariableNames',...
    {'freq', 'L0_good', 'L2_good','L4_good','L0_poor','L2_poor','L4_poor',...
    'L0_good_corr', 'L2_good_corr','L4_good_corr','L0_poor_corr','L2_poor_corr','L4_poor_corr'});

%% plots

cols = [204, 0, 0; ...
        0, 204, 0]/255;

figure
strtPlot = [1, 4, 7];
str_load = {'load0', 'load2', 'load4'};


big_poor = cat(3,HR_poor_L0, HR_poor_L2, HR_poor_L4);
big_good = cat(3,HR_good_L0, HR_good_L2, HR_good_L4);


for iPlot = 1:3
    
    StructLoad(iPlot).good = apply_detrend(squeeze(big_good(:,:,iPlot)), params);
    StructLoad(iPlot).poor = apply_detrend(squeeze(big_poor(:,:,iPlot)), params);
   
        
    ax = subplot(3,3,[strtPlot(iPlot),strtPlot(iPlot)+1]);
    lineProps.col = {cols(1,:)};
    lineProps.transparent = .5;
    
    vectAVG_poor = mean(StructLoad(iPlot).poor,2);
    vectAVG_good = mean(StructLoad(iPlot).good,2);
    
%     [StructLoad(iPlot).FIT.poor.freq, StructLoad(iPlot).FIT.poor.gof] = ...
%         createFitSin(params.time_bins, vectAVG_poor);
%     
%     [StructLoad(iPlot).FIT.good.fitresult, StructLoad(iPlot).FIT.good.gof] = ...
%         createFitSin(params.time_bins, vectAVG_good);
    
    mseb(params.time_bins', mean(StructLoad(iPlot).poor,2)',...
        std(StructLoad(iPlot).poor,[],2)'/sqrt(size(StructLoad(iPlot).poor,2)),...
        lineProps,1)
    hold on
    lineProps.col = {cols(2,:)};

    mseb(params.time_bins', mean(StructLoad(iPlot).good,2)',...
        std(StructLoad(iPlot).good,[],2)'/sqrt(size(StructLoad(iPlot).good,2)),...
        lineProps,1)
    
    plot(params.time_bins, zeros(numel(params.time_bins),1), 'k', 'LineWidth',.5)
    if iPlot ==1
        title({'detrended HR', str_load{iPlot}})
    
    else
        title(str_load{iPlot})
        
    end
    ylabel('detrended HR')
     
    if iPlot<3
        set(ax(1),'XTickLabel','')
        set(ax(1),'XColor',get(gca,'Color'))
        set(ax(1),'box','off')
        
         
    else
        
        xlabel('time (s)')
        
    end
    
    xlim(minmax(params.time_bins'))
 
    
    
      
end

% load 0
ax = subplot(3,3,3);

lineProps.col = {cols(2,:)};
mseb(current_spectra_good.freqs, AVG_spctr_good_L0',...
        stdERR_spctr_good_L0',...
        lineProps,1)
hold on

lineProps.col = {cols(1,:)};
mseb(current_spectra_good.freqs, AVG_spctr_poor_L0',...
        stdERR_spctr_poor_L0',...
        lineProps,1)
hold on


plot(spctr_good_L0.freqs, upper_good_L0, '--', 'Color', cols(2,:), 'LineWidth', 2); hold on
plot(spctr_good_L0.freqs, upper_poor_L0, '--','Color', cols(1,:), 'LineWidth', 2); hold on
%xlabel('freq (Hz)')
ylabel('fft amplitude')
xlim([2 10])
set(ax(1),'XTickLabel','')
set(ax(1),'XColor',get(gca,'Color'))
set(ax(1),'box','off')
        


title({'spectra and','permutations'})

% load 2

ax = subplot(3,3,6);


lineProps.col = {cols(2,:)};
mseb(current_spectra_good.freqs, AVG_spctr_good_L2',...
        stdERR_spctr_good_L2',...
        lineProps,1)
hold on

lineProps.col = {cols(1,:)};
mseb(current_spectra_good.freqs, AVG_spctr_poor_L2',...
        stdERR_spctr_poor_L2',...
        lineProps,1)
hold on

plot(spctr_good_L2.freqs, upper_good_L2, '--', 'Color', cols(2,:), 'LineWidth', 2); hold on
plot(spctr_good_L2.freqs, upper_poor_L2, '--','Color', cols(1,:), 'LineWidth', 2); hold on
%xlabel('freq (Hz)')
ylabel('fft amplitude')
% legend('good performers', 'good performers 95° percentile', 'poor performers',...
%     'poor performers 95° percentile','Location','eastoutside')
xlim([2 10])
%title('load 2')
set(ax(1),'XTickLabel','')
set(ax(1),'XColor',get(gca,'Color'))
set(ax(1),'box','off')
        

% load 4

subplot(3,3,9)


lineProps.col = {cols(2,:)};
mseb(current_spectra_good.freqs, AVG_spctr_good_L4',...
        stdERR_spctr_good_L4',...
        lineProps,1)
hold on

lineProps.col = {cols(1,:)};
mseb(current_spectra_good.freqs, AVG_spctr_poor_L4',...
        stdERR_spctr_poor_L4',...
        lineProps,1)
hold on

plot(spctr_good_L0.freqs, upper_good_L4, '--', 'Color', cols(2,:), 'LineWidth', 2); hold on
plot(spctr_good_L0.freqs, upper_poor_L4, '--','Color', cols(1,:), 'LineWidth', 2); hold on
xlabel('freq (Hz)')
ylabel('fft amplitude')
xlim([2 10])
%title('load 4')



%% try fitting approach


    
[GOF_perm.good.mat, GOF_perm.poor.mat] = deal(nan(n_iter, 3));
    
% obtain adjR2 out of permutation
for iPerm = 1:n_iter
    
    for iLoad = 1:3
        
        good_avg = mean(squeeze(perm_AND_det.good(:,:,iLoad, iPerm)),2);
        poor_avg = mean(squeeze(perm_AND_det.poor(:,:,iLoad, iPerm)),2);
        
        
        [~,GOF_perm.good.mat(iPerm, iLoad)] = createFitSin(params.time_bins,...
            good_avg); 
        [~,GOF_perm.poor.mat(iPerm, iLoad)] = createFitSin(params.time_bins,...
            poor_avg);
        
    end
    
end

GOF_perm.good.thresh95 = prctile(GOF_perm.good.mat, 95);
GOF_perm.poor.thresh95 = prctile(GOF_perm.poor.mat, 95);

figure
for iPlot = 1:3
    
    StructLoad(iPlot).good = apply_detrend(squeeze(big_good(:,:,iPlot)), params);
    StructLoad(iPlot).poor = apply_detrend(squeeze(big_poor(:,:,iPlot)), params);
    
        
    ax = subplot(3,1,iPlot);
    lineProps.col = {cols(1,:)};
    lineProps.transparent = .5;
    
    vectAVG_poor = mean(StructLoad(iPlot).poor,2);
    vectAVG_good = mean(StructLoad(iPlot).good,2);
    
    [StructLoad(iPlot).FIT.poor.freq, StructLoad(iPlot).FIT.poor.gof] = ...
        createFitSin(params.time_bins, vectAVG_poor);
    
    [StructLoad(iPlot).FIT.good.fitresult, StructLoad(iPlot).FIT.good.gof] = ...
        createFitSin(params.time_bins, vectAVG_good);
    
    mseb(params.time_bins', mean(StructLoad(iPlot).poor,2)',...
        std(StructLoad(iPlot).poor,[],2)'/sqrt(size(StructLoad(iPlot).poor,2)),...
        lineProps,1)
    hold on
    lineProps.col = {cols(2,:)};

    mseb(params.time_bins', mean(StructLoad(iPlot).good,2)',...
        std(StructLoad(iPlot).good,[],2)'/sqrt(size(StructLoad(iPlot).good,2)),...
        lineProps,1)
    
    
    
    
      
end

%% compare sin fit with 95 perc of perms
[mat_fitted_adjR2, best_freqfit, perc_95_perms] = deal(nan(3, 2));


for icond = 1:3
    
    mat_fitted_adjR2(icond, 1) = StructLoad(icond).FIT.poor.gof;
    best_freqfit(icond, 1) = StructLoad(icond).FIT.poor.freq;
    perc_95_perms(icond, 1) = GOF_perm.poor.thresh95(icond);

    mat_fitted_adjR2(icond, 2) = StructLoad(icond).FIT.good.gof;
    best_freqfit(icond, 2) = StructLoad(icond).FIT.good.fitresult;
    perc_95_perms(icond, 2) = GOF_perm.good.thresh95(icond);

end
