%% BEHAVIOURAL DATA ANALYSIS
% started: 21-04-2018
%

clear all
close all
clc

% logical flags
want_label = false;
want_all_trls = true;
raw_merged_already = true;


%% select raw data and merge it, or directly load the merged datafile
if ~raw_merged_already


    %% DEFINE FOLDERS RAWDATA

    % select datafile
    [fileVar,filePath] = uigetfile('*','choose the datafile','MultiSelect','on');
    if filePath==0, error('None selected!'); end



    %% GATHER ALL INFO INTO MATRIX FOR ALL SUBJECTS

    n_subject = numel(fileVar);


    mat_data = nan(7,1440,n_subject);

    for iSubj = 1:n_subject

        load([filePath fileVar{iSubj}])

        %% CREATE SINGLE MATRIX gathering all conditions, response and so on
        % 1 row: load condition
        % 2 row: mem equality condition
        % 3 row: subj mem response
        % 4 row: deltaT
        % 5 row: timestamp delta
        % 6 row: flash_presence
        % 7 row: subj flash response

        mat_data(1,:,iSubj) = main.condLoad;
        mat_data(2,:,iSubj) = main.eqVSdiff;
        mat_data(3,:,iSubj) = main.resp_MEM;
        mat_data(4,:,iSubj) = main.vect_deltaT;
        mat_data(5,:,iSubj) = main.timestamp(:,5)';
        mat_data(6,:,iSubj) = main.flashVScatch;
        mat_data(7,:,iSubj) = main.resp_FLASH;


    end

else
    
    load('../RAW_data_collapsed.mat')
    n_subject = size(mat_data, 3);
    
end

%% START ANALYSIS

[ACC_mat, HR_mat ]= deal(nan(3,16,n_subject));

% define SDT functions (either parametric and non-parametric)
make_dPrime = @(HR,FAR) norminv(HR,0,1)-norminv(FAR,0,1);

make_criter = @(HR,FAR) -.5*(norminv(HR,0,1)+norminv(FAR,0,1));

make_Aprime = @(HR,FAR) .5+(sign(HR-FAR)*((HR-FAR)^2+abs(HR-FAR))...
    /(4*max([HR, FAR])-4*HR*FAR)); % generalized formula from todorov 1999

make_Bsecond = @(HR,FAR) sign(HR-FAR)*(HR*(1-HR)-FAR*(1-FAR))...
    /(HR*(1-HR) + FAR*(1-FAR)); % as before, non parametric measure of decision bias ased on todorov 1999

SDT_mat = nan(2,2,3,n_subject);   % HR  |  MR   || load || subject
                                  % FAR |  CRR  ||

[dPrime_mat, criterion_mat, Aprime_mat, Bsecond_mat] = ...
    deal(nan(n_subject,3));

[mem_accuracy, cowanK_mat] = deal(nan(n_subject,2));
                                  
for iSubj = 1:n_subject

    curr_mat = mat_data(:,:,iSubj);
        
    % warning --nan + get delta values 
    delta_vals = unique(curr_mat(4,:));
    if any(isnan(delta_vals))
        error(['NaN in delta vals for subj ' num2str(iSubj) '. not good'])
    end
    
    load_vals = unique(curr_mat(1,:));
    
    
    %% start loop for conditions
    
    % find accurate memory trials
    memAcc_Lmask = curr_mat(2,:)==curr_mat(3,:);
    memAcc_Lmask(curr_mat(1,:)==0)=1; % set to correct trials with 0 load
    
    wrongTrlsMEM(:,iSubj) = abs(memAcc_Lmask-1);
    
    loopLoad = 0;
    for iLoad = load_vals
        
        loopLoad = loopLoad+1;
        
        if loopLoad>1
            
            innerLoop = loopLoad-1;
            mem_accuracy(iSubj,innerLoop) = mean(memAcc_Lmask(curr_mat(1,:)==iLoad));
            
            % Cowan's K computation
            % K = N*(HR+CR-1)
            mem_EQUAL_vect = (curr_mat(3,:)==1); % subj response "equal"
            mem_EQUAL = mean(mem_EQUAL_vect((curr_mat(1,:)==iLoad)&curr_mat(2,:)));
            
            mem_DIFF_vect = (curr_mat(3,:)==0);
            mem_DIFF = mean(mem_DIFF_vect((curr_mat(1,:)==iLoad)&(curr_mat(2,:)==0)));
            
            cowanK_mat(iSubj, innerLoop) = iLoad*(mem_DIFF+mem_EQUAL-1);
            
        end
    
        % logical vector for load condition
        load_Lmask = curr_mat(1,:)==iLoad;
        
        % logical vector for load and memory accuracy
        loadMemAcc_Lmask = load_Lmask&memAcc_Lmask;
        
        % select whole dataset regardless of memory accuracy
        if want_all_trls
            
            loadMemAcc_Lmask = load_Lmask;
            
        end
        
        %% SDT computation X condition       
        reduced_mat_SDT = curr_mat([6 7], loadMemAcc_Lmask);
        
        % logical masks for flash present or absent, respectively
        vect_flash = reduced_mat_SDT(1,:)==1;
        vect_catch = reduced_mat_SDT(1,:)==0;
        
        % Hit Rate
        hits_vect = reduced_mat_SDT(2,vect_flash);
        HR = nanmean(hits_vect);
        SDT_mat(1,1,loopLoad,iSubj) = HR;
        
        % Miss Rate
        miss_vect = abs(hits_vect-1);
        MR = nanmean(miss_vect);
        SDT_mat(1,2,loopLoad,iSubj) = MR;
        
        % False Alarm Rate
        falseAl_vect = reduced_mat_SDT(2,vect_catch);
        if all(falseAl_vect==0) % FAR correction
            falseAl_vect(1) = 1;
        end
        FAR = nanmean(falseAl_vect);
        SDT_mat(2,1,loopLoad,iSubj) = FAR;

        % Correct Rejections
        corrRej_vect = abs(falseAl_vect-1);
        if all(corrRej_vect==1) % CRR correction
            corrRej_vect(1) = 0;
        end
        CRR = nanmean(corrRej_vect);
        SDT_mat(2,2, loopLoad, iSubj) = CRR;
        
        % dPrime
        dPrime_mat(iSubj, loopLoad) = make_dPrime(HR, FAR);
        
        % criterion
        criterion_mat(iSubj, loopLoad) = make_criter(HR, FAR);
        
        % Aprime 
        Aprime_mat(iSubj, loopLoad) = make_Aprime(HR, FAR);
        
        % Bsecond
        Bsecond_mat(iSubj, loopLoad) = make_Bsecond(HR, FAR);
        
        
        %% START LOOP FOR DELTA T
        
        loopDelta = 0;
        for iDelta = delta_vals
            
            loopDelta = loopDelta+1;
            
            % find current delta t in data
            deltaT_Lmask = curr_mat(4,:)==iDelta; disp(sum(deltaT_Lmask));

            % mask for current load, delta t in right memory trials
            delta_load_mem_Lmask = deltaT_Lmask&loadMemAcc_Lmask;
            
            % mask for HR (only flash present)
            flash_present = curr_mat(6,:)==1;
            HR_lMask = delta_load_mem_Lmask&flash_present;
            
            % put HR into matrix, for each load, each delta, each
            % participant
            curr_HR = nanmean(curr_mat(7,HR_lMask));
            HR_mat(loopLoad, loopDelta, iSubj) = curr_HR;
            
            % try accuracy as well
            current_delta_acc_mat = curr_mat([6 7], delta_load_mem_Lmask);
            current_delta_acc_vect = current_delta_acc_mat(1,:)==...
                current_delta_acc_mat(2,:);
            curr_acc = nanmean(current_delta_acc_vect);
            ACC_mat(loopLoad,loopDelta,iSubj) = curr_acc;
            
            
            
            
            
            
            
        end

        
    end
        
  
    
end

%% ACROSS SUBJECTS AVERAGE SDT
ga_SDT = squeeze(mean(SDT_mat,4));
stdERR_SDT = squeeze(std(SDT_mat,0, 4))/sqrt(n_subject);

ga_dPrime = mean(dPrime_mat);
stdERR_dPrime = std(dPrime_mat)/sqrt(n_subject);

ga_criterion = mean(criterion_mat);
stdERR_crit = std(criterion_mat)/sqrt(n_subject);

ga_Aprime = mean(Aprime_mat);
stdERR_Aprime = std(Aprime_mat)/sqrt(n_subject);

ga_Bsecond = mean(Bsecond_mat);
stdERR_Bsecond = std(Bsecond_mat)/sqrt(n_subject);


%% PLOT SDT RESULTS
str_titles = {'Hit Rate (HR)', 'False Alarm (FA)', 'd''', 'criterion'};
str_loads = {'load0','load2','load4'};
cols = [0 204 204;
        127 0 255;
        255 51 153]/255;


figure

%% HR
curr_plot = 1;
subplot(2,2,curr_plot)

ss_hr = squeeze(SDT_mat(1,1,:,:))';
mask_plot = repmat(1:3,n_subject,1);


for iScatter =1:3
    
    scatter(mask_plot(:,iScatter)+randn(n_subject,1)/20, ss_hr(:,iScatter),[],cols(iScatter,:), 'filled')
    
    if want_label
        text(mask_plot(:,iScatter), ss_hr(:,iScatter), fileVar)
    end
    
    hold on 
    
end

[HR_anova, HR_stdERR, HR_MC] = rm1W_ANOVA_adapted(ss_hr, str_loads, 0,0,'')


errorbar(1:3, squeeze(ga_SDT(1,1,:)), HR_stdERR, 'k',...
    'Linewidth', 3)
hold on

set(gca, 'XTickLabel',str_loads, 'XTick',1:3)
title(str_titles{curr_plot});

%% FA
curr_plot = 2;
subplot(2,2,curr_plot)

ss_fa = squeeze(SDT_mat(2,1,:,:))';


for iScatter =1:3
    
    scatter(mask_plot(:,iScatter)+randn(n_subject,1)/20, ss_fa(:,iScatter),[], cols(iScatter,:),'filled')
    
    if want_label
        text(mask_plot(:,iScatter), ss_fa(:,iScatter), fileVar)       
    end
    
    hold on 
    
end

[FA_anova, FA_stdERR, FA_MC] = rm1W_ANOVA_adapted(ss_fa, str_loads, 0,0,'')

errorbar(1:3, squeeze(ga_SDT(2,1,:)), FA_stdERR, 'k',...
    'Linewidth', 3)
hold on

set(gca, 'XTickLabel',str_loads, 'XTick',1:3)
title(str_titles{curr_plot});

%% dPrime
curr_plot = 3;
subplot(2,2,curr_plot)


for iScatter =1:3
    
    scatter(mask_plot(:,iScatter)+randn(n_subject,1)/20, dPrime_mat(:,iScatter),[], cols(iScatter,:),'filled')
    
    if want_label
        text(mask_plot(:,iScatter), dPrime_mat(:,iScatter), fileVar)       
    end
    
    hold on 
    
end

[dPrime_anova, dPrime_stdERR, dPrime_MC] = rm1W_ANOVA_adapted(dPrime_mat, str_loads, 0,0,'')

errorbar(1:3, ga_dPrime, stdERR_dPrime, 'k',...
    'Linewidth', 3)
hold on

set(gca, 'XTickLabel',str_loads, 'XTick',1:3)
title(str_titles{curr_plot});


%% criterion
curr_plot = 4;
subplot(2,2,curr_plot)


for iScatter =1:3
    
    scatter(mask_plot(:,iScatter)+randn(n_subject,1)/20, criterion_mat(:,iScatter),[], cols(iScatter,:),'filled')
    
    if want_label
        text(mask_plot(:,iScatter), criterion_mat(:,iScatter), fileVar)       
    end
    
    hold on 
    
end

[criterion_ANOVA, stdERR_crit, crit_MC] = rm1W_ANOVA_adapted(criterion_mat, str_loads, 0,0,'')


errorbar(1:3, ga_criterion, stdERR_crit, 'k',...
    'Linewidth', 3)
hold on

set(gca, 'XTickLabel',str_loads, 'XTick',1:3)
title(str_titles{curr_plot});



%% check SDT in good vs poor

thresh_GvsP = median(cowanK_mat(:,2));
idx_good = find(cowanK_mat(:,2)>=thresh_GvsP);
idx_poor = find(cowanK_mat(:,2)<thresh_GvsP);

%% plot decrease VWM perf
figure
errorbar(1:2, mean(mem_accuracy), std(mem_accuracy)/sqrt(n_subject),...
    'k','LineWidth',2)

hold on
scatter(ones(numel(idx_good),1)+randn(numel(idx_good),1)/20,...
    mem_accuracy(idx_good,1), [],[0 204/255 0],'filled')

scatter(ones(numel(idx_poor),1)+randn(numel(idx_poor),1)/20,...
    mem_accuracy(idx_poor,1), [], [204/255 0 0],'filled')

scatter(2*ones(numel(idx_good),1)+randn(numel(idx_good),1)/20,...
    mem_accuracy(idx_good,2), [],[0 204/255 0],'filled')

scatter(2*ones(numel(idx_poor),1)+randn(numel(idx_poor),1)/20,...
    mem_accuracy(idx_poor,2), [], [204/255 0 0],'filled')



set(gca, 'XTickLabel',{'load2', 'load4'}, 'XTick',1:2)
title('VWM accuracy')
ylim([.6 1])



HR_good = ss_hr(idx_good,:);
HR_poor = ss_hr(idx_poor,:);

FA_good = ss_fa(idx_good,:);
FA_poor = ss_fa(idx_poor,:);

criterion_good = criterion_mat(idx_good,:);
criterion_poor = criterion_mat(idx_poor,:);

dPrime_good = dPrime_mat(idx_good,:);
dPrime_poor = dPrime_mat(idx_poor,:);

SDT_GvsP = [];
%HR
SDT_GvsP.HR.good.avg =  mean(HR_good);
SDT_GvsP.HR.good.sderr = std(HR_good)/sqrt(length(HR_good));
SDT_GvsP.HR.poor.avg = mean(HR_poor);
SDT_GvsP.HR.poor.sderr = std(HR_poor)/sqrt(length(HR_poor));
%FA
SDT_GvsP.FA.good.avg =  mean(FA_good);
SDT_GvsP.FA.good.sderr = std(FA_good)/sqrt(length(FA_good));
SDT_GvsP.FA.poor.avg = mean(FA_poor);
SDT_GvsP.FA.poor.sderr = std(FA_poor)/sqrt(length(FA_poor));
% dprime
SDT_GvsP.dPrime.good.avg =  mean(dPrime_good);
SDT_GvsP.dPrime.good.sderr = std(dPrime_good)/sqrt(length(dPrime_good));
SDT_GvsP.dPrime.poor.avg = mean(dPrime_poor);
SDT_GvsP.dPrime.poor.sderr = std(dPrime_poor)/sqrt(length(dPrime_poor));
% criterion
SDT_GvsP.criterion.good.avg =  mean(criterion_good);
SDT_GvsP.criterion.good.sderr = std(criterion_good)/sqrt(length(criterion_good));
SDT_GvsP.criterion.poor.avg = mean(criterion_poor);
SDT_GvsP.criterion.poor.sderr = std(criterion_poor)/sqrt(length(criterion_poor));

nFields = fieldnames(SDT_GvsP)';

c = 0;
figure
for iField = nFields
    
    c=c+1;
    subplot(2,2,c)
    hold on
    
    errorbar(1:3, SDT_GvsP.(iField{1}).good.avg, ...
        SDT_GvsP.(iField{1}).good.sderr, 'Color', [0 204/255 0],'LineWidth', 2)
    
    errorbar(1.15:3.15, SDT_GvsP.(iField{1}).poor.avg, ...
        SDT_GvsP.(iField{1}).poor.sderr, 'Color', [204/255 0 0],'LineWidth', 2)
    
    set(gca, 'XTickLabel',str_loads, 'XTick',1:3)

    title(iField)
    
    if c == 4
        legend('good performers','poor performers')
        
    end
    
end


%% save output

save('../ACC_allsubj.mat','ACC_mat')
save('../CowansK.mat', 'cowanK_mat')
save('../dPrime.mat', 'dPrime_mat')

if want_all_trls
    
    save('../HR_allsubj_allTrls.mat','HR_mat')
    
else
    
    save('../HR_allsubj.mat','HR_mat')

end
