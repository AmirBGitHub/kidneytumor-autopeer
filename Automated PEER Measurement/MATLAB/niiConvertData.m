clear all; close all; clc

for i = 1:100
    img_CT = permute(readNPY(['./data/CT_img_subj',num2str(i),'.npy']), [3 2 1]);
    niftiwrite(img_CT, ['./data/niiData/CT_img_nii_subj',num2str(i),'.nii']);
    img_KT = permute(readNPY(['./data/KT_img_subj',num2str(i),'.npy']), [3 2 1]);
    niftiwrite(img_KT, ['./data/niiData/KT_img_nii_subj',num2str(i),'.nii']);
end
