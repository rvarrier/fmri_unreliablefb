clearvars;
dbstop if error;


addpath('G:\Backup\ODT\algorithms\spm12\spm12');

data_loc = 2;
device = 2;
if device == 1
    rootfile_data = 'E:\ODT';
else
    rootfile_data = 'G:\Backup\ODT';
end
addpath(fullfile(rootfile_data,'algorithms','cvmanova-3')); % https://github.com/allefeld/cvmanova


data_folder_name = pwd;
cd(data_folder_name);

[subj_nm, num_subj, num_bl, start_bl, end_bl, nvol,nvol_prepost, nvol_int,nvol_floc,start_blnum_perrun,coreg_bw_runs] = all_subj_names();

subj_list = 6:35;

manovaType = 2; % manovaType = input('1. searchlight or 2. roi-based?: ');

four_regressors = 1; % sep. regressors for LL, LR, RL, RR
if four_regressors
    spm_fileloc = 'fourReg';
else
    spm_fileloc = 'twoReg';
end
colin_extn = '_colin';

roi_str = {'V1';'V123';'LIP';'MFG'};
mask_folder_name = {'Visual_hOc1';'Visual_hOc123'; roi_str{3}; roi_str{4}};


if manovaType == 1
    slRadius = 3; % 21mm diameter
end
setnum = [1 3 4];
subj_ind = 0;

for subj = subj_list % subj_list = 6 to 35
    subj_ind = subj_ind + 1;
    if subj <= 13;        num_runs = 6; % 2 runs each removed from each test phase in 8 subjects
    else;        num_runs = 8;     end
    Cs = {[1 -1]'
        [0 0 1 -1]'}; % F-contrasts. Corresponding GLM has first four regressors as CCW135,CW135, CCW45 and CW45 for each run.
    
    for set_ind = 1:2 % pre/post-intervention
        
        if manovaType == 1 % searchlight analysis
            spm_filename = strcat('spm_run',num2str(setnum(set_ind)), '_norm_sl'); % normalised space
            dirName = fullfile(data_folder_name, [num2str(subj) '.' subj_nm{subj}], 'fmri', spm_fileloc, spm_filename);
            cvManovaSearchlight(dirName, slRadius, Cs);
            
        else % ROI analysis
            spm_filename = strcat('spm_run',num2str(setnum(set_ind)), '_colin'); % colin space
            thr = 0.5; % probability of voxels to belong to V1 from SPM anatomy toolbox
            spm_floc = ['spm_run4' '_colin'];
            mask_folder = fullfile(data_folder_name, [num2str(subj) '.' subj_nm{subj}], 'fmri', mask_folder_name{roi});
            cd(mask_folder);
            mask_file = ls(['*_anatMask' num2str(thr) '.nii']);
            cd(data_folder_name);
            %                 if size(mask_file,1) == 0 % anat file needs to be created
            loc_with_maindata = 0;
            [~, roi_flname] = extract_anat_visual_voxels(subj, thr, data_folder_name, spm_fileloc, spm_floc,...
                four_regressors);
            
            dirName = fullfile(data_folder_name, [num2str(subj) '.' subj_nm{subj}], 'fmri', spm_fileloc, spm_filename);
            
            [D, vox] = cvManovaRegion(dirName, roi_flname, Cs); % D = 2*1 matrix (2 contrasts: CCW135 vs. CW135 and CCW45 vs. CW45)
            D_all{subj_ind,set_ind} = D; % each cell element is a 2*1 matrix
            vox_all(subj_ind,set_ind) = vox;
        end
    end
end



if manovaType == 2
    thr = 0.5; % prob of V1
    thr_str = num2str(thr);
    thr_extn = ['_thr_' thr_str(1) '_' thr_str(3)];
    
    save(fullfile(data_folder_name, 'cvMANOVA_results', ['cvMANOVA' thr_extn '_' roi_str{roi} '_colin' '.mat']), 'D_all', 'vox_all');
    
end