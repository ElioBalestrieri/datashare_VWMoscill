function [outTableANOVA, within_StdErr, outTableMC] = rm1W_ANOVA_adapted(data,VariableNames,saveOut,graph,titlegraph)
%% [outTable] = rm1W_ANOVA(data,VariableNames)
%  
% data: row--> subj; col-->variables of interest
% VariableNames must be a cell containing strings with dim== n° columns of data
% saveOut--> if saveOut set to 1, save out to txt file; 
% graph--> if==1 show barplot
% titlegraph --> string with graph title
%
% wrote by Nicholas Menghi
% adapted by Elio Balestrieri and Nicholas Menghi [25-Sep-2017]
% added multiple comparisons (eb) sometime in 2018 summer
% added eta square computation (eb) [16-Feb-2019]

% Set the value
temp = size(data);
n = temp(1); %How many subjects per condition
k = temp(2); % How many condition
N = temp(1)*temp(2); % Total subjects (even if repeated)

%% one-way Anova - dependent samples analysis (Repeated measure)

% conditional mean
cndl_means = mean(data);
%grandERRb = std(data)/sqrt(n);
    
% grand mean
grand_mean = mean(cndl_means);

%% SS within
% error variation of individual score around each group mean
% not a product of manipulation
%temp = data-repmat(cndl_means,n,1);

temp = [];
for i = 1:k
    temp(:,i) = data(:,i) - cndl_means(i);
end


ss_within = sum(sum(temp.^2));

%% SS between
% Variance among group
ss_between = sum(((cndl_means-grand_mean).^2)*n); % why times n???

%% SS Participant
% compute sum of square per participant mean
participant_means = mean(data,2);
ss_participant = sum((participant_means-grand_mean).^2*k);

%% SS error
% Error of the model
ss_error = ss_within - ss_participant;

%% ss total
ss_total = sum(sum((data-mean(mean(data))).^2));

eta_squared = ss_between/ss_total;

%degree of freedom
df_participant = n-1;
df_between = k-1;
df_error = (n-1)*(k-1);
df_total = N-1;

% ms
ms_participant = ss_participant ./ df_participant;
ms_between = ss_between ./ df_between;
ms_error = ss_error ./ df_error;

F = ms_between ./ ms_error;
CumulativeProbability = fcdf(F, df_between, df_error); % Cumulative probability till f Value
p =  1- fcdf(F, df_between, df_error);
% 1/0 cause is HIGHLY significant

%% calculate error bars

for zz = 1:size(data,2)
    
    meanDataErrorBars(:,zz) = data(:,zz)-mean(data,2);
end

meanDataErrorBars = meanDataErrorBars+mean(mean(data));
within_StdErr = std(meanDataErrorBars)/sqrt(size(meanDataErrorBars,1));


%% create figure
if graph==1
    % Create new figure
    %h = figure;
        
    % Plot bar graph with mean per column represented
    bar(mean(data),'FaceColor',[0.1 0.40 .90])
    hold on
    errorbar(cndl_means,within_StdErr,'.k')
    set(gca,'fontsize',14)
    set(gca, 'XTickLabel',VariableNames, 'XTick',1:numel(VariableNames),'FontSize',20)
    %set(gca,'TitleFontSizeMultiplier',2)
    %set(gca,'LabelFontSizeMultiplier',2)
    %ylabel('FontSize',20)
    title(titlegraph, 'FontSize',24)
    ylim([0 max(cndl_means)+max(cndl_means)/3])
    %legend(['ANOVA p = ' num2str(p)])
end


%% plot
% figure;
% x = -0.05:0.01:6;
% y = fpdf(x, df_between, df_error); %Probability density function
% % gives back the probability for each value of x
% plot(x, y, 'b');

%% table
outTableANOVA=table(df_participant,df_total, ms_participant, F, CumulativeProbability,p, eta_squared);
%fprintf('\n\n')
%disp(outTable)

%% multiple comparison

combs_mc = nchoosek(1:numel(VariableNames),2);
n_comp = size(combs_mc,1);

outTableMC = [];

for iMC = 1:n_comp
    
    outTableMC(iMC).var1 = VariableNames{combs_mc(iMC,1)};
    outTableMC(iMC).var2 = VariableNames{combs_mc(iMC,2)};
    
    smpl1 = data(:,combs_mc(iMC,1));
    smpl2 = data(:,combs_mc(iMC,2));
    
    [~,P,CI,STAT] = ttest(smpl1, smpl2);
    
    outTableMC(iMC).CI = CI;
    outTableMC(iMC).t = STAT.tstat;
    outTableMC(iMC).df = STAT.df;
    outTableMC(iMC).pVal_uncorrected = P;
    outTableMC(iMC).pVal_bonferroni = P*n_comp;
    outTableMC(iMC).Cohens_d = (mean(smpl1)-mean(smpl2))/...
        sqrt((var(smpl1)+var(smpl2))/2);    
    
end
    

if isempty(outTableMC)==0
    outTableMC = struct2table(outTableMC);
end

if saveOut==1
    write(outTableANOVA,[ titlegraph '_' datestr(now)  'RM1W_ANOVA.csv']);
end

end

