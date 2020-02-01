% code to analyse eyetracking data accuracy and multivariate analyses.

clearvars;
dbstop if error;
dbstop if warning;

data_folder_name =  'C:\Work\data_backup_eye\eye_data_GitHub';
cd(data_folder_name);
load('data_eyetracking.mat');
margin = 68; % fixation window in pxls
margin_deg = 1.3; % fixation window in degrees
max_x = 1024; % screen res x
max_y = 768; % screen res y

%% sanity check
row = 0;
for fb = 1:2
    for subj = 1:size(stim_data_per_fbtype_phase_subj,2)
        row = row + 1;
        num_elements = [];
        for run = [1 3] % pre-/post-int test phases
            if run == 1;                run_ind = 1;
            else; run_ind = 2; end
            num_elements(run_ind) = size(stim_data_per_fbtype_phase_subj{fb,subj,run},1);
            dat = stim_data_per_fbtype_phase_subj{fb,subj,run};
            if size(dat,1) > 744 % if data file has at least 30% data
                pc_null(fb,subj,run) = 100.*numel(find(isnan(dat(:,6))))/size(dat(:,6),1); % compute invalid data pts (in %)  in file 
            else
                pc_null(fb,subj,run) = 100; % excluded rightaway
            end
            
            r_bl = [];
            
            if pc_null(fb,subj,run) ~= 100
                x = dat(:,6); y = dat(:,7); % x and y coord in pxls (pre-processed and baseline-corrected)
                bl_vec = dat(:,4); % block numbers (1-8)
                invalid_x = (x >= max_x) | (x <= 0) | isnan(x);
                invalid_y = (y >= max_y) | (y <= 0) | isnan(y);
                invalid_xy = invalid_x | invalid_y;
                
                % remove NaNs
                x = x(~invalid_xy); 
                if ~isempty(find(isnan(x)))
                    error('nanx');
                end
                y = y(~invalid_xy);
                if ~isempty(find(isnan(y)))
                    error('nany');
                end
                bl_vec = bl_vec(~invalid_xy);
                
                r = sqrt((x - 0.5*max_x).^2 + (y - 0.5*max_y).^2); % distance from the origin
                dat1 = dat(~invalid_xy,:);
                
                for bl = 1:8
                    rows_bl = find(bl_vec == bl);
                    r_bl = r(bl_vec == bl);
                    if ~isempty(find(isnan(r_bl)))
                        error('nan');
                    end
                    valid_r(bl) = numel(find(r_bl <= margin));
                    total_r(bl) = numel(r_bl);
                end
                pc_fix(row,run_ind) = 100*sum(valid_r)/sum(total_r); % %fixation
            
            else % too many missing data pts
                pc_fix(row,run_ind) = nan;
            end
        end
        pc_fix(row,3) = fb;
        
        if  ~isempty(find(pc_null(fb, subj,[1 3]) >= 70)) % if at least one test phase has missing data
            subj_incl(row) = 0;
            pc_fix(row,1:2) = nan(2,1);
        
        else subj_incl(row) = 1;
        end
        
    end
end

pc_fix_sel = pc_fix(subj_incl ==1,:);
rows1 = find(pc_fix_sel(:,3) == 1);
rows2 = find(pc_fix_sel(:,3) == 2);
mean_pc_fix = [mean(pc_fix_sel(rows1,1:2)); mean(pc_fix_sel(rows2,1:2))];
std_pc_fix = [std(pc_fix_sel(rows1,1:2)); std(pc_fix_sel(rows2,1:2))];
se_pc_fix = [std(pc_fix_sel(rows1,1:2))/sqrt(numel(rows1));...
    std(pc_fix_sel(rows2,1:2))/sqrt(numel(rows2))];


%% multivariate

subj_incl = [];
addpath(genpath('C:\Work\libsvm-3.24\'));
load('data_across_blocks_and_subjects_incl_no_resp.mat'); % including no response trials 

row = 0;
for fb = 1:2
    subj_list_perfbtype = subj_list{fb};
    for subj = 1:size(subj_list_perfbtype,1)
        row = row + 1; 
        num_elements = [];
        
        if  ~isempty(find(squeeze(pc_null(fb, subj,[1 3])) >= 70)) % not enough valid data pts in file
            subj_incl(row) = 0;
        else
            subj_incl(row) = 1;
        end
        
        subj_ind_behdata = subj_num_beh{fb,subj};
            
        for run = [1 3]
            if run == 1;                run_ind = 1;
            else; run_ind = 2; end
            dat = stim_data_per_fbtype_phase_subj{fb,subj,run}; % cols: 1.fbtype,2.subjind(1-10rel. fb,11-22unrel. fb)...
            % 3.run(pre/post), 4.blno (1 to 8), 5.trialnr.(2-32), 6.x(inclNaNs), 7.y(inclNaNs), 8.pupilsize_x, 9.pupilsize_y
            
            %             flName = fullfile('C:\Work\ODT\ODT\data\Analysis',subj_folder,'behavioural',[subj_ind '_TeBl' ...
%                 num2str(run_ind) '.mat']),'data_testing_runs';
%             beh = load(flName);

            stim_type = [];% to store block, trialnum and stimtype from 32*8 trials
            for bl = 1:8
                col = (run_ind-1)*16+bl;
                dat_per_block = Aalldata_incl_no_resp{subj_ind_behdata,col};
                stim_type =[stim_type;dat_per_block(:,1:3)];
            end
            
            if run == 1
                thresh(row) = thr(subj_ind_behdata);
                subj_id(row) = subj_ind_behdata;
            end
            
            if size(dat,1) > 744
                blno = dat(:,4);
                stimno = dat(:,5);
                it = 0;
                x = dat(:,6); y =dat(:,7); % pre-processed bsln corrected x and y pixel values (incl NaNs)
                ctr_x = 0.5*max_x; ctr_y = 0.5*max_y; % pixel values of the centre of the screen
                
                dist_from_ctr = [x-ctr_x y-ctr_y];
                
                % MVPA
                for bl = 1:8
                    for stim = 2:32
                        it = it + 1;
                        rows_perstim_perbl = find(dat(:,4) == bl & dat(:,5) == stim);
                        for i = 1:length(rows_perstim_perbl)
                            xp = dat(rows_perstim_perbl(i),8); % pupil_size x (baseline-corrected and roundness (%diff between x and y) > 0.9)
                            yp = dat(rows_perstim_perbl(i),9); % pupil_size y (baseline-corrected and roundness (%diff between x and y) > 0.9)
                            if (xp <= 0) && (yp <= 0) % consttriction
                                xyp(i) = -(sqrt(double((xp.^2) + (yp.^2))));
                            elseif (xp >= 0) && (yp >= 0) % dilation
                                xyp(i) = sqrt(double((xp.^2) + (yp.^2)));
                            elseif ((xp > 0) && (yp < 0)) || ((xp < 0) && (yp > 0))
                                error('x and y diff signs');
                            end
                        end
                        
                        sel_row_stim = find((stim_type(:,1) == bl) & (stim_type(:,2)==stim));
%                         if ~isempty(sel_row) 
                            feature_matrix = [nanmean(dat(rows_perstim_perbl,6)) nanstd(dat(rows_perstim_perbl,6)) ...
                                nanmean(dat(rows_perstim_perbl,7))  nanstd(dat(rows_perstim_perbl,7))...
                                nanmean(xyp) nanstd(xyp)]; % nanmean(velocity_dat_deg(rows))];
                            pattern_per_subj_phase(it,1:9) = [feature_matrix stim_type(sel_row_stim,3) stim bl];
                            pattern_per_subj_phase(pattern_per_subj_phase == 0) = NaN;
%                         end
                    end % stim
                end % bl
                
                for r = 1: size(pattern_per_subj_phase,1) % approx 32*8 rows
                    %                             if ~isempty(find(isnan(pattern_per_subj_phase(r,1:10))))
                    nan_perrow = find(isnan(pattern_per_subj_phase(r,:)));
                    if numel(nan_perrow) == 0
                        nanrow(r) = 0;
                    elseif numel(nan_perrow) > 2 % if more than 2 ourt of 6 features are missing, remove the row
                        nanrow(r) = 1;
                    else % if <=2  NaNs out of 6 features, impute it if there are >80% elements  (i.e. 80% of 31 trials within the block)...
                        % for that feature
                        for c = 1:numel(nan_perrow) % e.g. 2 nans
                            nan_ind = nan_perrow(c);
                            if numel(find(isnan(pattern_per_subj_phase(:,nan_ind)))) ...
                                    < 0.8*size(pattern_per_subj_phase,1) % <80% nan in that column
                                pattern_per_subj_phase(r,nan_ind) = nanmedian(pattern_per_subj_phase(:,nan_ind)); % imputation
                            end
                            %                                     disp('wait');
                        end
                        nan_perrow = find(isnan(pattern_per_subj_phase(r,:)));
                        if numel(nan_perrow) > 0 % if imputation did not help remove all NaNs
                            nanrow(r) = 1;
                        end
                    end % nans in a row
                end % row 1 to all rows
                pattern_per_subj_phase = pattern_per_subj_phase(~nanrow,:); % remove rows with missing features even after imputaiton
                
                if ~isempty(find(isnan(pattern_per_subj_phase)))
                    error('nans');
                end
                
                all_bl = zeros(8,1);
                % to make this comparable to the fMRI data, only 6 blocks are selected for subj 1-8 
                if subj_num_beh{fb,subj} <= 8
                    bl_list = [1:4 6 8];
                else
                    bl_list = 1:8;
                end
                all_bl(bl_list) = 1;
                
                incl_row = [];
                for r = 1:size(pattern_per_subj_phase,1)
                    if all_bl(pattern_per_subj_phase(r,9)) == 1
                        incl_row(r) = 1;
                    else
                        incl_row(r) = 0;
                    end
                end
                mat = pattern_per_subj_phase(incl_row==1,:); % excludes bl5 and 7 already
                
                for bl_ind = 1: length(bl_list) % 1 to 6 or 1 to 8, n-1 cross-validation
                    disp(['test bl:' num2str(bl_ind)]);
                    mat_train = mat(mat(:,9) ~= bl_list(bl_ind),:); % ideally 31*7 =217 trials
                    mat_test = mat(mat(:,9) == bl_list(bl_ind),:); % ideally 31 trials
                    if ((size(mat_train,1) > 217)||(size(mat_test,1)>31))
                        error('check data');
                    end
                    
                    train_1 = (mat_train(mat_train(:,7) == 1,1:6)); % ideally around 54 trials, stimtype1
                    train_2 = (mat_train(mat_train(:,7) == 2,1:6)); % ideally around 54 trials, stimtype2
                    train_3 = (mat_train(mat_train(:,7) == 3,1:6)); % ideally around 54 trials, stimtype3
                    train_4 = (mat_train(mat_train(:,7) == 4,1:6)); % ideally around 54 trials, stimtype4
                    test_1 = (mat_test(mat_test(:,7) == 1,1:6)); % ideally around 7.75 trials
                    test_2 = (mat_test(mat_test(:,7) == 2,1:6)); % ideally around 7.75 trials
                    test_3 = (mat_test(mat_test(:,7) == 3,1:6)); % ideally around 7.75 trials
                    test_4 = (mat_test(mat_test(:,7) == 4,1:6)); % ideally around 7.75 trials
                    
                    train_set_pair1 = [train_1;train_2];
                    test_set_pair1 = [test_1;test_2];
                    labels_train_pair1 = [ones(size(train_1,1),1);2*ones(size(train_2,1),1)];
                    labels_test_pair1 = [ones(size(test_1,1),1);2*ones(size(test_2,1),1)];
                    if (size(train_set_pair1,1) > 10) && (size(test_set_pair1,1) > 1)  % too few rows
                        mdl_libsvm_pair1 = svmtrain(labels_train_pair1,train_set_pair1,['-q -t 0 -c 1' ]);
                        [~,acc_libsvm_pair1,~] = svmpredict(labels_test_pair1, test_set_pair1, mdl_libsvm_pair1,'-q');
%                         DA12 = 100.*numel(find(labels_test-lbls_num))/size(lbls,1);
                        DA12 = acc_libsvm_pair1(1);
                    else
                        DA12 = NaN;
                    end
                    
                    train_set_pair2 = [train_3;train_4];
                    test_set_pair2 = [test_3;test_4];
                    labels_train_pair2 = [3*ones(size(train_3,1),1);4*ones(size(train_4,1),1)];
                    labels_test_pair2 = [3*ones(size(test_3,1),1);4*ones(size(test_4,1),1)];
                    if (size(train_set_pair2,1) > 10) && (size(test_set_pair2,1) > 1)
                        mdl_libsvm_pair2 = svmtrain(labels_train_pair2,train_set_pair2,['-q -t 0 -c 1' ]);
                        [~,acc_libsvm_pair2,~] = svmpredict(labels_test_pair2, test_set_pair2, mdl_libsvm_pair2,'-q');
%                         DA34 = 100.*numel(find(labels_test-lbls_num))/size(lbls,1);
                        DA34 = acc_libsvm_pair2(1);
                    else
                        DA34 = NaN;
                    end
                    DA(bl_ind) = mean([DA12 DA34]);
                    wt(bl_ind) = numel(find(~isnan([DA12 DA34])))/2;
                    
                end
             mean_DA(row,run_ind) =  nansum(DA.*wt)/sum(wt); %weighted sum depending on how much data was available for each iteration
                
            else % not enough data pts in the file
                pc_null(fb,subj,run) = 100;
                mean_DA(row,run_ind) = NaN;
            end
        end
        
        if ~isempty(find(isnan(mean_DA(row,:)))) || ...
                (~isempty(find(mean_DA(row,:) == -1))) % missing test phases
            mean_DA (row,run_ind) = NaN;
        end
        mean_DA(row,3) = fb;
        mean_DA(row,4) = subj_ind_behdata;
    end
end


%% MeanDA
figure; hold on; xlim([0 3]);
set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[1 1 8 6]);

par = {'*-b'; 'o-r'};
subj_incl_col = subj_incl(:);
for fb = 1:2
    rows_perstim_perbl = find((mean_DA(:,3) == fb) & (subj_incl_col == 1));
    if fb == 1;
        dat_plot1_DA = mean_DA(rows_perstim_perbl,1:2);
        thr11 = thresh(rows_perstim_perbl);
        subj_id1 = subj_id(rows_perstim_perbl);
        diff_dat_plot1 = dat_plot1_DA(:,2) - dat_plot1_DA(:,1);
    else;
        thr22 = thresh(rows_perstim_perbl);
        subj_id2 = subj_id(rows_perstim_perbl);
        dat_plot2_DA = mean_DA(rows_perstim_perbl,1:2);
        diff_dat_plot2 = dat_plot2_DA(:,2) - dat_plot2_DA(:,1);
    end
end
addpath(genpath('C:\Work\klabhub-bayesFactor-V2.0-0-gba51949'));
[bf10,p,CI,stats] = bf.ttest2(diff_dat_plot1,diff_dat_plot2);
M = [dat_plot1_DA ones(size(dat_plot1_DA,1),1) thr11' subj_id1';...
    dat_plot2_DA 2*ones(size(dat_plot2_DA,1),1) thr22' subj_id2'];
csvwrite(fullfile(data_folder_name,'Decoding_accuracy_backup.csv'),M);
mean_dat_plot_DA = [mean(dat_plot1_DA); mean(dat_plot2_DA)];
ste_dat_plot_DA = [std(dat_plot1_DA)./sqrt(size(dat_plot1_DA,1)); ...
    std(dat_plot2_DA)./sqrt(size(dat_plot2_DA,1))];
[~, p_pre_fb1] = ttest(dat_plot1_DA(:,1)-50);
[~, p_post_fb1] = ttest(dat_plot1_DA(:,2)-50);
[~, p_pre_fb2] = ttest(dat_plot2_DA(:,1)-50);
[~, p_post_fb2] = ttest(dat_plot2_DA(:,2)-50);

bar(1, mean_dat_plot_DA(1,1),'FaceColor', 'none', 'EdgeColor',[0 0 0]);
bar(2, mean_dat_plot_DA(1,2),'FaceColor', [0.5 0.5 0.5], 'EdgeColor',[0 0 0]);
bar(4, mean_dat_plot_DA(2,1),'FaceColor', 'none', 'EdgeColor',[0 0 0]);
bar(5, mean_dat_plot_DA(2,2),'FaceColor', [0.5 0.5 0.5], 'EdgeColor',[0 0 0]);

plot([1 2], dat_plot1_DA,'Color',[0.3 0.3 0.3],'LineWidth', 0.5);
plot([4 5], dat_plot2_DA,'Color',[0.3 0.3 0.3],'LineWidth', 0.5);
plot(0:6, repmat(50,1,7), '--k');
legend('Pre-intervention','Post-intervention', 'Location', 'northeast');
%         xlim([0.8 2.2]);
xlim([0.5 6]);
ylim([30 70]);
ylabel('Mean DA');
set(gca,  'FontSize', 7, 'XTick',[1.5 4.5], 'XTickLabel',{'Reliable fb', 'Unreliable fb'});
savefig(fullfile(data_folder_name,'MeanDA'));
saveas(gcf,fullfile(data_folder_name,'MeanDA'),'jpg');
print(gcf,'-dpng','-r300',fullfile(data_folder_name,'MeanDA'));