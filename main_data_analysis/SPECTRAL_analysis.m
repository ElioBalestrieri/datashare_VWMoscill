function [outANOVAS, structT] = SPECTRAL_analysis(L0, L2, L4, x_freq)

% reduced power spectra
reduced_powerspectraHR = cat(3, L0, L2, L4);
n_subj = size(reduced_powerspectraHR,2);
str_load = {'load 0', 'load 2', 'load 4'};

%% compute ANOVAs
[TTESTS_0vs4,TTESTS_0vs2,TTESTS_2vs4,ANOVAS, outANOVAS ]= deal(nan(numel(x_freq), 5));
 

for iFreq = 1:numel(x_freq)
    
    [outTable,~,~] = rm1W_ANOVA_adapted(squeeze(reduced_powerspectraHR(iFreq,:,:)),...
        [],0,0,[]);
    
    outANOVAS(iFreq,1) = x_freq(iFreq);
    outANOVAS(iFreq,2) = .1;
    outANOVAS(iFreq,3) = outTable.F;
    outANOVAS(iFreq,4) = outTable.p;

    
    if outTable.p <.05

        ANOVAS(iFreq,1) = x_freq(iFreq);
        ANOVAS(iFreq,2) = .01;
        ANOVAS(iFreq,3) = outTable.F;
        ANOVAS(iFreq,4) = outTable.p;

    end
    
    % 0vs4
    [TTESTS_0vs4(iFreq, 4), TTESTS_0vs4(iFreq, 2),~, stat_t]...
        = ttest(squeeze(reduced_powerspectraHR(iFreq,:,1))',...
        squeeze(reduced_powerspectraHR(iFreq,:,3))');
    
    TTESTS_0vs4(iFreq,3) = stat_t.tstat;
    TTESTS_0vs4(iFreq,1) = x_freq(iFreq);
    TTESTS_0vs4(iFreq,5) = TTESTS_0vs4(iFreq,2)*numel(x_freq);
    
    % 0vs2
    [TTESTS_0vs2(iFreq, 4), TTESTS_0vs2(iFreq, 2),~, stat_t]...
        = ttest(squeeze(reduced_powerspectraHR(iFreq,:,1))',...
        squeeze(reduced_powerspectraHR(iFreq,:,2))');
    
    TTESTS_0vs2(iFreq,3) = stat_t.tstat;
    TTESTS_0vs2(iFreq,1) = x_freq(iFreq);
    TTESTS_0vs2(iFreq,5) = TTESTS_0vs2(iFreq,2)*numel(x_freq);
    
    % 2vs4
    [TTESTS_2vs4(iFreq, 4), TTESTS_2vs4(iFreq, 2),~, stat_t]...
        = ttest(squeeze(reduced_powerspectraHR(iFreq,:,2))',...
        squeeze(reduced_powerspectraHR(iFreq,:,3))');
    
    TTESTS_2vs4(iFreq,3) = stat_t.tstat;
    TTESTS_2vs4(iFreq,1) = x_freq(iFreq);
    TTESTS_2vs4(iFreq,5) = TTESTS_2vs4(iFreq,2)*numel(x_freq);


end

structT = [];
structT.L0vsL2 = TTESTS_0vs2;
structT.L2vsL4 = TTESTS_2vs4;
structT.L0vsL4 = TTESTS_0vs4;


% AVG HR spectra
AVG_pwrspectra = mean(reduced_powerspectraHR,2);
SDER_pwrspectra = std(reduced_powerspectraHR, 0, 2)/sqrt(n_subj);

figure
colors = get(gca, 'ColorOrder');
for iPlot = 1:3
    lineProps.col = {colors(iPlot,:)};
    mseb(x_freq, AVG_pwrspectra(:,iPlot)', SDER_pwrspectra(:,iPlot)',lineProps,1)
    
end
hold on
scatter(ANOVAS(:,1), ANOVAS(:,2), 40, [0 1 0], 'filled')

legend(str_load)
title(['Grand Avg spectra HR --' num2str(n_subj) ' subj--'])
xlim([min(x_freq) max(x_freq)])

 
end