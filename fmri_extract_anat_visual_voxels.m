function [anatMask, roi_flname] = extract_anat_visual_voxels(subj,thr,data_folder_name,spm_fileloc,spm_floc,...
    four_regressors)
%% source of data and job files
thr_extn = num2str(thr);
mask_folder_name = {'Visual_hOc1';'Visual_hOc123'};
mask_file_name = {'wVisual_hOc1';['wVisual_hOc123_' thr_extn(1) '_' thr_extn(3)]};

job_folder_name =  data_folder_name;
cd(job_folder_name);
[subj_nm, num_subj, num_bl, start_bl, end_bl, nvol,nvol_prepost, nvol_int,nvol_floc,start_blnum_perrun,coreg_bw_runs] = all_subj_names();

%%
folder_name = fullfile(data_folder_name,[num2str(subj) '.' subj_nm{subj}],'fmri', mask_folder_name{1}, [mask_file_name{1} '.nii']);

anatMask_all = spm_read_vols(spm_vol(folder_name));
maskName = fullfile(job_folder_name,[num2str(subj) '.' subj_nm{subj}], 'fmri',spm_fileloc,spm_floc,'mask.nii');
mask = spm_read_vols(spm_vol(maskName));

if four_regressors
    folder_name_visloc = fullfile(data_folder_name,[num2str(subj) '.' subj_nm{subj}],'fmri',spm_fileloc,spm_floc,'spmT_0005.nii');
else
    folder_name_visloc = fullfile(data_folder_name,[num2str(subj) '.' subj_nm{subj}],'fmri',spm_fileloc,spm_floc,'spmT_0003.nii');
end

visualMask = spm_read_vols(spm_vol(folder_name_visloc));
visualMask_thr = zeros(size(visualMask));

load(fullfile(data_folder_name,[num2str(subj) '.' subj_nm{subj}],'fmri',spm_fileloc,spm_floc,'SPM.mat'));
df = SPM.xX.erdf;
thr_visual = tinv(0.95,df); 

visualMask_thr(find(visualMask> thr_visual)) = 1; % voxels common to LL, Lr,RL,RR, i.e., task-relevant voxels

anatMask_all = anatMask_all .* mask .* visualMask_thr; % anatMask_all = anatomical, mask = mask from floc giving valid voxels, 
                                          % visualMask_thr = voxels >95%
                                          % thr (p < 0.05)
anatMask_all(find(isnan(anatMask_all))) = 0;
% end

anatMask = zeros(size(anatMask_all));
anatMask(find(anatMask_all > thr)) = 1; % anatMask = 1 if it's (1) above threshold, (2) a valid voxel from mask.nii of floc
    % and (3) valid floc (all stim >0, p 0.001 unc.)


cd(fullfile(data_folder_name,[num2str(subj) '.' subj_nm{subj}],'fmri',mask_folder_name{1}));
sourceFile = folder_name;
destFile = strcat(sourceFile(1:end-4),'_anatMask',num2str(thr), '.nii'); % this name was later changed to anatMask_Snum_....
copyfile(sourceFile, destFile);
V = spm_vol(destFile);
spm_write_vol(V, anatMask);
roi_flname = destFile;
cd(job_folder_name);


