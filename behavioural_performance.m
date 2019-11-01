clear all;
dbstop if error;

device = 2;
data_loc = 2; % data_loc = input('Data in 1. S-drive or 2. ext HD: ');
if data_loc == 1
    rootfile_data = 'S:\AG\AG-Sterzer\Varrier\ODT';
else
    if device == 1
        rootfile_data = 'E:\Backup\ODT';
    else
        rootfile_data = 'G:\Backup\ODT';
    end
end

data_folder_name = fullfile(pwd,'beh_results');
cd(data_folder_name);
load('data_across_blocks_and_subjects');

num_subj = size(Aalldata,1);
missing_values = zeros(num_subj,1);
for subj = subj_list % subj_list = 6:35
    for run_ind = 1:24
        data_per_bl = Aalldata{subj,run_ind};
        start_trial = 1;
        data_per_bl = data_per_bl(start_trial:end,:); % select data
        
        %% define DVs
        stim = data_per_bl(:,3); stim_lr = mod(stim,2); stim_lr(stim_lr==0) = 2;
        resp_corrected = data_per_bl(:,4);  resp_lr = resp_corrected;
        if resp_lr > 2; error('wrong resp'); end
        num_trials = size(data_per_bl,1);
        num_trials_all(subj,run_ind) = num_trials;
        
        %DV1
        correct_ind = find(stim_lr == resp_lr); wrong_ind = find(stim_lr ~= resp_lr);
        num_correct_all(subj,run_ind) =  length(correct_ind);
        correct_ind_l = find(stim_lr == resp_lr & resp_lr == 1); num_trials_l = length(find(stim_lr == 1));
        correct_ind_r = find(stim_lr == resp_lr & resp_lr == 2); num_trials_r = length(find(stim_lr == 2));
        percperf(subj,run_ind) = 100.*( length(correct_ind)/num_trials );
        
    end
    % missed resp
    for set_ = 1:3 % baseline test phase, intervention phase, post-intervention test phase
        runs_ = (set_-1)*8+1:set_*8;
        missed_resp_perset(subj,set_) = mean(num_missed_resp(subj,runs_));
    end
end

num_trials_all = num_trials_all(subj_list,:);
num_correct_all = num_correct_all(subj_list,:);
percperf = percperf(subj_list,:); % 30*24
Afbtype = Afbtype(subj_list);% 1*30
Afbtype = Afbtype';
thr = thr(subj_list);% 1*30
missed_resp_perset = missed_resp_perset(subj_list,:);
subj_sel_missedtrials = subj_sel_missedtrials(subj_list);
tbl_DV1 = table(percperf,Afbtype,thr);
save('tbl_DV1.mat','tbl_DV1');
if device ==1
    writetable(tbl_DV1,fullfile(data_folder_name ,'DV1_percperf.csv'));
end
meanperf = 100.*[sum(num_correct_all(:,1:8),2)./sum(num_trials_all(:,1:8),2) ...
    sum(num_correct_all(:,9:16),2)./sum(num_trials_all(:,9:16),2)...
    sum(num_correct_all(:,17:24),2)./sum(num_trials_all(:,17:24),2)];

meanperf_diff = meanperf(:,3) - meanperf(:,1);
for fb = 1:2
    meanperf_diff_perfbtype{fb} = meanperf_diff(find(Afbtype == fb));
end

subj_sel = ones(30,1);
subj_noref = ones(30,1);
% subj_noref([8 13 27]) = 0;
subj_sel = subj_noref;
disp(['num subj: ' num2str(numel(find(subj_sel)))]);
cont = input('continue? 1. yes, 0. no: ');
if cont ~= 1
    error('stop, change num subj');
end

tbl_DV1_mean = table(meanperf,Afbtype,thr,missed_resp_perset,subj_sel_missedtrials, subj_noref);
DV1_mean = [meanperf Afbtype thr missed_resp_perset subj_sel_missedtrials subj_noref];
save('tbl_DV1_mean.mat','tbl_DV1_mean');
writetable(tbl_DV1_mean,fullfile(data_folder_name ,'DV1_percperf_mean.csv'));

%% to compare behavioural and fmri data, some runs need to be excluded
for subj = 1:30
    if subj <=8
        runs = [1:4 6 8];
    elseif subj == 18
        runs = [1:2 4:8];
    else runs = 1:8;
    end
    
    avg_fmri(subj,1) = 100.*sum(num_correct_all(subj,runs),2)./sum(num_trials_all(subj,runs),2);
    avg_fmri(subj,2) = 100.*sum(num_correct_all(subj,16+runs),2)./sum(num_trials_all(subj,16+runs),2);
    
    %     avg_fmri_old(subj,1) = mean(percperf(subj,runs));
    %     avg_fmri_old(subj,2) = mean(percperf(subj,16+runs));
    
end
tbl_DV1_mean_fmri = table(avg_fmri,Afbtype,thr,missed_resp_perset,subj_sel_missedtrials);
save('tbl_DV1_mean_fmri.mat','tbl_DV1_mean_fmri');
writetable(tbl_DV1_mean_fmri,fullfile(data_folder_name ,'DV1_percperf_mean_fmri.csv'));


%                 fun_percperf(DV1);
%         handle = figure('units','normalized','outerposition',[0 0 1 1]);
%         numrows = 7;
%         numcol = 5;
%         if numrows*numcol < num_subj
%             error('change dimensions of subplot');
%         end
%         for subj = 1:30
%             missing_values(subj) = numel(find((percperf(subj,:)==0)| isnan(percperf(subj,:))));
%             curr_rownum = ceil(subj/numcol);
%             curr_colnum = mod(subj,numcol); curr_colnum(curr_colnum == 0) = numcol;
%             left = (curr_colnum-1)*(1/numcol) + 0.1*(1/numcol); width = 0.8*(1/numcol);
%             bottom = (numrows-curr_rownum).*(1/numrows) + 0.1*(1/numrows); height = 0.8*(1/numrows);
%             axes('Position',[left bottom width height]);
%             %             subplot(numrows,numcol,subj);
%             hold on; axis([0 25 40 100]);
%             title(strcat('subj ',num2str(subj)));
%             plot(1:8,percperf(subj,1:8),'.k','MarkerSize',6,'LineWidth', 1.5);
%             plot(1:8, repmat(mean(percperf(subj,1:8)),1,8),'--k');
%             plot(17:24,percperf(subj,17:24),'.k','MarkerSize',6,'LineWidth', 1.5);
%             plot(17:24, repmat(nanmean(percperf(subj,17:24)),1,8),'--k');
%             if Afbtype(subj) == 1
%                 plot(9:16,percperf(subj,9:16),'*b','MarkerSize',6,'LineWidth', 1.5);
%                 plot(9:16, repmat(nanmean(percperf(subj,9:16)),1,8),'--b');
%             else
%                 plot(9:16,percperf(subj,9:16),'or','MarkerSize',6,'LineWidth', 1.5);
%                 plot(9:16, repmat(nanmean(percperf(subj,9:16)),1,8),'--r');
%             end
%         end
%         savefig(fullfile(data_folder_name,'Fig1'));
%         saveas(gcf,fullfile(data_folder_name,'Fig1'),'jpg');


handle = figure('units','normalized','outerposition',[0.25 0.25 0.3 0.5]);
hold on; set(gca,'FontSize',11);
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[1 1 7 7]);
axis([0 4 40 100]);
colour = {'*-b';'o-r'};
for fb = 1:2
    %             ind = find(Afbtype==fb);
    ind = find(Afbtype==fb & subj_sel == 1);
    div = sqrt(numel(ind));
    errorbar(mean(meanperf(ind,[1 3])), std(meanperf(ind,[1 3]))/div, colour{fb},'MarkerSize',6,'LineWidth', 1.5);
    
end

legend('Correct fb', 'Unreliable fb','Location','northwest');
%         legend('boxoff');
ylabel('% correct');
set(gca,'XTick',1:2,'XTickLabel',{'Pre';'Post';})
xlim([0.8 2.2]);
ylim([70 90]);
savefig(fullfile(data_folder_name,'Fig2'));
saveas(gcf,fullfile(data_folder_name,['Fig2_n' num2str(numel(find(subj_sel)))]),'jpg');
print(gcf,'-dpng','-r300',fullfile(data_folder_name,['Fig2_n' num2str(numel(find(subj_sel)))]));



%% Main figure plot bars withiout errorbars
figure; hold on;
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[1 1 10 7]);

for fb = 1:2
    %             ind = find(Afbtype==fb);
    ind = find(Afbtype==fb & subj_sel == 1);
    div = sqrt(numel(ind));
    mean_perf_per_fbtype{fb} = meanperf(ind,[1 3]);
end
mean_mean_perf = [nanmean(mean_perf_per_fbtype{1}); nanmean(mean_perf_per_fbtype{2})];
ste_mean_perf = [nanstd(mean_perf_per_fbtype{1})./sqrt(size(mean_perf_per_fbtype{1},1));...
    nanstd(mean_perf_per_fbtype{2})./sqrt(size(mean_perf_per_fbtype{2},1))];

bar(1, mean_mean_perf(1,1),'FaceColor', 'none', 'EdgeColor',[0 0 0]);
bar(2, mean_mean_perf(1,2),'FaceColor', [0.5 0.5 0.5], 'EdgeColor',[0 0 0]);
bar(4, mean_mean_perf(2,1),'FaceColor', 'none', 'EdgeColor',[0 0 0]);
bar(5, mean_mean_perf(2,2),'FaceColor', [0.5 0.5 0.5], 'EdgeColor',[0 0 0]);

errorbar([1 2], mean_mean_perf(1,:), ste_mean_perf(1,:),'.k','LineWidth', 1);
errorbar([4 5], mean_mean_perf(2,:), ste_mean_perf(2,:),'.k','LineWidth', 1);

plot([1 2], mean_perf_per_fbtype{1},'Color',[0.3 0.3 0.3],'LineWidth', 0.5);
plot([4 5], mean_perf_per_fbtype{2},'Color',[0.3 0.3 0.3],'LineWidth', 0.5);

legend('Pre-intervention','Post-intervention', 'Location', 'northeast');
set(gca, 'FontSize', 8, 'XTick', 1:2, 'XTickLabel', {'Pre', 'Post'});
%         xlim([0.8 2.2]);
xlim([0.5 6]);
ylim([40 110]);
ylabel('Percent correct');
set(gca, 'XTick',[1.5 4.5], 'XTickLabel',{'Reliable fb', 'Unreliable fb'}, 'YTick',40:10:100);
savefig(fullfile(data_folder_name,'Fig2_subjwise'));
saveas(gcf,fullfile(data_folder_name,['Fig2_subjwise_n' num2str(numel(find(subj_sel)))]),'jpg');
print(gcf,'-dpng','-r300',fullfile(data_folder_name,['Fig2_subjwise_n' num2str(numel(find(subj_sel)))]));
%%
[delta_fb1.h, delta_fb1.p, ~, delta_fb1.stats] = ttest(meanperf_diff_perfbtype{1});
mean_meanperf_diff_fb1 = mean(meanperf_diff_perfbtype{1});
se_meanperf_diff_fb1 = std(meanperf_diff_perfbtype{1})/sqrt(numel(meanperf_diff_perfbtype{1}));
[delta_fb2.h, delta_fb2.p, ~, delta_fb2.stats] = ttest(meanperf_diff_perfbtype{2});
mean_meanperf_diff_fb2 = mean(meanperf_diff_perfbtype{2});
se_meanperf_diff_fb2 = std(meanperf_diff_perfbtype{2})/sqrt(numel(meanperf_diff_perfbtype{2}));

std1 = std(meanperf_diff_perfbtype{1});
std2 = std(meanperf_diff_perfbtype{2});
std_pooled = sqrt( 0.5*( (std1^2) + (std2^2) ) );
cohens_D = (mean_meanperf_diff_fb2 - mean_meanperf_diff_fb1)/std_pooled;
N = 30;
f = ((N-3)./(N-2.25)).* sqrt((N-2)/N);
cohens_D_corrected = cohens_D.*f;

disp('debug');

%% plot bars with errorbars
figure; hold on;
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[1 1 9 7]);
bar(1, mean_mean_perf(1,1),'FaceColor', 'none', 'EdgeColor',[0 0 0]);
bar(2, mean_mean_perf(1,2),'FaceColor', [0.5 0.5 0.5], 'EdgeColor',[0 0 0]);
bar(4, mean_mean_perf(2,1),'FaceColor', 'none', 'EdgeColor',[0 0 0]);
bar(5, mean_mean_perf(2,2),'FaceColor', [0.5 0.5 0.5], 'EdgeColor',[0 0 0]);
for fb = 1:2
    for i = 1:15
        s = ind(i);
        perc_perf_per_fbtype{i,fb} = percperf(s,[1:8 17:24]);
    end
end

for fb = 1:2
    for s = 1:15
        dat = perc_perf_per_fbtype{s,fb};
        den = sqrt(length(dat)/2);
        xloc = (fb-1)*3+1:(fb-1)*3+2;
        errorbar(xloc,[nanmean(dat(:,1:8)) nanmean(dat(:,9:16))],[nanstd(dat(:,1:8)) nanstd(dat(:,9:16))]./den,'k')
    end
end

legend('Pre-intervention','Post-intervention', 'Location', 'northeast');
set(gca, 'FontSize', 8, 'XTick', 1:2, 'XTickLabel', {'Pre', 'Post'});
%         xlim([0.8 2.2]);
xlim([0.5 5.5]);
ylim([20 120]);
ylabel('Percent correct');
set(gca, 'XTick',[1.5 4.5], 'XTickLabel',{'Reliable fb', 'Unreliable fb'});
savefig(fullfile(data_folder_name,'Fig2_subjwise_errorbars'));
saveas(gcf,fullfile(data_folder_name,['Fig2_subjwise_errorbars_n' num2str(numel(find(subj_sel)))]),'jpg');
print(gcf,'-dpng','-r300',fullfile(data_folder_name,['Fig2_subjwise_errorbars_n' num2str(numel(find(subj_sel)))]));