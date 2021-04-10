%% permute and sinusoid fit

clearvars
close all
clc

want_all_trls = false;

if want_all_trls
    
    load('../HR_allsubj_allTrls.mat')
    load('../HR_allsubj_L02_alltrls.mat')

else
    
    load('../HR_allsubj.mat')
    load('../HR_allsubj_L02.mat')

end

% reproducibility
rng(0)


%
n_subj = size(HR_mat,3);
n_tp = size(HR_mat,2);


% merge together load 02 and load 4
expnd_L02 = reshape(HR_L02, 1, n_tp, n_subj);
HR_mat_lowhigh = cat(1, expnd_L02, HR_mat(3, :, :));


n_iter = 100;

str_load = {'low load', 'high load'};

detrended_HR = nan(size(HR_mat_lowhigh));
detrended_shuffled_HR = nan(n_tp, 2, n_subj, n_iter);

% permute from raw, and then detrend inside the big loop
HR_mat_perm = permute_from_raw_lowhighload(want_all_trls, n_iter);

params = [];
params.detrend_flag = 2;
params.time_bins = (.15:.04:.750)';

for iLoad = 1:2
    
    for iSubj = 1:n_subj
        
        [detrended_HR(iLoad,:,iSubj),~,~] = ...
            do_detrend(squeeze(HR_mat_lowhigh(iLoad,:,iSubj))', 1 , params);
        
        for iIter = 1:n_iter
            
            real_vect = HR_mat(iLoad,:,iSubj);
            rsmpld_vect = HR_mat_perm(iLoad, :, iSubj, iIter);
            [detrended_shuffled_HR(:,iLoad,iSubj,iIter),~,~] = ...
                do_detrend(rsmpld_vect', 1, params);
            
            
       end
        
    end
    
    
    
end



AVG_HR_det_real = mean(detrended_HR,3)';

SDERR_HR_det_real = std(detrended_HR,[],3)'/sqrt(n_subj);

AVG_HR_PermAndDet = squeeze(mean(detrended_shuffled_HR, 3));

%% FIT SINUSOID

h_fit_full = fittype( 'sin1' ); 

x = (.150:.040:.750);
          
lowerBands = [0;
              0;
              0];
upperBands = [Inf;
              Inf;
              Inf];
          
          
permuted_adj_R2 = nan(n_iter, 2);
          
real_adj_R2 = nan(1,2);

figure

cell_load_fits = cell(1,3);

colors = [0 204 204;
         255 51 153]/255;

plotCount = [1, 4];

for iLoad = 1:2

    [xData, yData] = prepareCurveData( x,(AVG_HR_det_real(:,iLoad)));
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Robust',...
        'Bisquare');
    opts.Display = 'Off';
    
    opts.Lower = [-Inf lowerBands(iLoad,1) -Inf]; % -Inf lowerBands(iLoad,2) -Inf];
    
    opts.Upper = [+Inf upperBands(iLoad,1) +Inf]; % +Inf upperBands(iLoad,2) +Inf];
    
    [cell_load_fits{iLoad}, statR, out] = fit(xData, yData , h_fit_full, opts);
    
    currPlot = plotCount(iLoad);
    
    ax = subplot(2,3,[currPlot, currPlot+1]); 
    hold on
    
    lineProps.col = {colors(iLoad,:)};
    lineProps.transparent = .5;
    %lineProps.width = .01;
    lineProps.style = '.';
    
    mseb(xData', yData', SDERR_HR_det_real(:,iLoad)', lineProps)
    
    hold on;
    scatter(xData,yData, 20, colors(iLoad,:), 'filled')
    hold on;
    plot(cell_load_fits{iLoad})
    xlim(minmax(x))
    
    swap = coeffvalues(cell_load_fits{iLoad});
    freqs_fitted = swap(2)/(2*pi); % 5]);
    str_title1 = ['Best frequency sinusoid fit: ' ...
        num2str(round(freqs_fitted(1),2)) ' Hz'];
    
    str_title2 = ['AdjR^2: ' num2str(round(statR.adjrsquare,2))];
    
    title({str_load{iLoad}, str_title1, str_title2});
    
    ylabel({'2nd order', 'detrended HR'})
    
    real_adj_R2(iLoad) = statR.adjrsquare;
    %real_adj_R2(iLoad) = statR.rsquare;
    
    if iLoad<2
        
        set(ax(1),'XTickLabel','')
        set(ax(1),'XColor',get(gca,'Color'))
        set(ax(1),'box','off')
        
    else 
        
        xlabel('time (s)')
    
        
    end
    
    
    for iIter = 1:n_iter
        
        [xData, yData] = prepareCurveData( x,...
            squeeze(AVG_HR_PermAndDet(:,iLoad, iIter))');
        
        [~, statP, ~] = fit(xData, yData , h_fit_full, opts);
        
        permuted_adj_R2(iIter, iLoad) = statP.adjrsquare;
        %permuted_adj_R2(iIter, iLoad) = statP.rsquare;
        
        if mod(iIter, 100)==0
                
                fprintf('%d permuted fittings \n', iIter)
                
        end
            
        
    end
    
    
end

%% plot adj R2 in relation to histograms

for iPlot = 1:2
    
    ax = subplot(2,3,iPlot*3);
    
    % take matrix of permutation
    h = histogram(permuted_adj_R2(:,iPlot),100, 'FaceColor', [.1 .1 .1],...
        'Normalization', 'probability');
    hold on
    lineUp = max(h.Values); 
    
    % compute empirical cdf on the permutations 
    [curr_cdf(:,1), curr_cdf(:,2)] = ecdf(permuted_adj_R2(:,iPlot));
    % reverse the vector and find the first value exceeding alpha (.05)
    curr_threshold = curr_cdf(find(abs(curr_cdf(:,1)-1)<.05,1),2);
    
    % plot 95° percentile (set preferred value on y axis to make the line
    % coherent with graph
    plot([curr_threshold curr_threshold], [0 lineUp], '--k', 'LineWidth',2)
    
    hold on
    
    % real value obtained
    curr_real_R2 = real_adj_R2(iPlot);
    
    % plot real value obtained
    plot([curr_real_R2, curr_real_R2], [0 lineUp],'Color', colors(iPlot,:),...
        'LineWidth', 2);
    
    % see where your p is
    pVal = 1-curr_cdf(find(curr_cdf(:,2)>=curr_real_R2,1)-1,1);
    
    % add label for p value
    text(curr_real_R2, lineUp-lineUp/10, ['p = ' num2str(pVal)])
    
    if iPlot<2
        
        set(ax(1),'XTickLabel','')
        set(ax(1),'XColor',get(gca,'Color'))
        set(ax(1),'box','off')
        
    else
        
        xlabel('Adj R^2')
        set(ax(1),'box','off')
        
    end

    
    if iPlot==1
        title('permuted AdjR^2')
    end
    ylim([0 lineUp])
    xlim([-.3 1])
    ylabel('probability')
    
end

    


