clear all; close all; clc

for i = 1:100
    img_CT = permute(readNPY(['./python image data/CT_img_subj',num2str(i),'.npy']), [3 2 1]);
    niftiwrite(img_CT, ['./python image data/niiData/',num2str(i),'_CT.nii']);
    img_KT = permute(readNPY(['./python image data/KT_img_subj',num2str(i),'.npy']), [3 2 1]);
    niftiwrite(img_KT, ['./python image data/niiData/',num2str(i),'_Label.nii']);
end
