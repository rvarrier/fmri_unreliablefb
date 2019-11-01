clear all;
dbstop if error;
addpath('T:\Dokumente\MATLAB\spm12\');
data_loc = 2;
device = 2;
if data_loc == 1;
    rootfile_data = 'S:\AG\AG-Sterzer\Varrier\ODT';
else;
    if device == 1
        rootfile_data = 'F:\ODT';
    else
        rootfile_data = 'H:\Backup\ODT';
    end
end
addpath(fullfile(rootfile_data,'data','Analysis'));
curr_dir = fullfile(rootfile_data,'data','Analysis','re_analysis_fmri_April2019');
suff = {'_pre';'_post'};
[subj_nm, num_subj, num_bl, start_bl, end_bl, nvol,nvol_prepost, nvol_int,nvol_floc,start_blnum_perrun,coreg_bw_runs] = all_subj_names();

for subj = [1:21 23:30]
    disp(['subj ' num2str(subj)]);
    for ph = [1 2]
        ROI = fullfile(curr_dir, ['S' num2str(subj) '_anatMask'], 'voxels_Mask.nii');
        V = spm_vol(ROI);
        mask = spm_read_vols(V);
        mask_ind = find(mask);
        
        con_loc = fullfile(curr_dir, ['S' num2str(subj) suff{ph}], 'con_0001.nii');
        V = spm_vol(con_loc);
        con = spm_read_vols(V);
        
        con_selected = con(mask_ind);
        mean_act(subj,ph) = nanmean(con_selected);
        numVox(subj,ph) = numel(find(~isnan(con_selected)));
%         disp('debug');
    end
    
    str = fullfile(rootfile_data,'data','Analysis', [num2str(subj+5) '.' subj_nm{subj+5}], 'behavioural', ...
        [subj_nm{subj+5} '_TeBl1.mat']);
    load(str);
    thr(subj,1) = sig_threshold;
    
end
mean_act(22,1:2) = nan(1,2);
numVox(22,1:2) = nan(1,2);
thr(22,1) = NaN;

fbtype = repmat([2;1],15,1);
fbtype(22) = NaN;
subj_sel = ones(30,1); subj_sel(22) = 0;
fbtype = fbtype(subj_sel==1);
mean_act = mean_act(find(subj_sel),:);
numVox = numVox(find(subj_sel),:);
thr = thr(find(subj_sel));
csvwrite(fullfile(curr_dir,'results_univ_act.csv'), [mean_act numVox fbtype thr]);
save(fullfile(curr_dir,'results_univ_act.mat'),'mean_act', 'numVox', 'fbtype', 'thr');
cd(curr_dir);

diff_mean_act =  mean_act(:,2)-mean_act(:,1);
[h,p] = ttest2(diff_mean_act(fbtype ==1), diff_mean_act(fbtype ==2))

mean_mean_act = [mean(mean_act(fbtype ==1,:)); mean(mean_act(fbtype ==2,:))];
se_mean_act = [std(mean_act(fbtype ==1,:))/sqrt(size(mean_act(fbtype ==1,:),1)); ...
    std(mean_act(fbtype ==2,:))/sqrt(size(mean_act(fbtype ==2,:),1))];

mean_diff_mean_act = [mean(diff_mean_act(fbtype ==1)), mean(diff_mean_act(fbtype ==2))];
se_diff_mean_act = [std(diff_mean_act(fbtype ==1))/sqrt(numel(find(fbtype==1))) ...
    std(diff_mean_act(fbtype ==2))/sqrt(numel(find(fbtype==2)))];

