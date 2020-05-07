# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Amir Baghdadi
"""
print(' ')
print('CNN for Onco vs. ChRCC Identification')
print('ATLAS Group')
print('Roswell Park Comprehensive Cancer Center')
print(' ')

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
#import cv2
from collections import Counter
import glob
import re
import os

CT_list_all = []
tumor_list_all = []
kidney_list_all = []
CT_img_all = []
KT_img_all = []
nsliceCT = []
nsliceKT = []

# read image data
subj = 1
for root, dirs, files in os.walk('main_images'):
    imageNum_list = []
    CT_list = []
    tumor_list = []
    kidney_list = []
    CT_img = []
    KT_img = []

    # main image
#    print('root: ', root)
#    print('dirs: ', dirs)
#    print('files: ', files)
    if not dirs:
        print('Patient #', subj)
        dir = ''.join([root,'\*[0-9].jpg'])
        for filename in glob.glob(dir):
#            print(filename)
            lst = re.findall('view([0-9]+)', filename) # get the numeric part of name
            imageNum_list.append(lst[0])
            lst = re.search('view(.*).jpg', filename) # get the whole name
            CT_list.append(lst[0])
            CT_img.append(mpimg.imread(filename))
    #        imgplot = plt.imshow(mpimg.imread(filename))
    #        plt.show()
        for i in range(len(imageNum_list)):    
            # tumor image
            filename_T = ''.join(['view',imageNum_list[i],'-T.jpg'])
            lst_T = re.search('view(.*).jpg', filename_T)

            # kidney image
            filename_K = ''.join(['view',imageNum_list[i],'-K.jpg'])
            lst_K = re.search('view(.*).jpg', filename_K)

            if filename_T in files and not filename_K in files:
#                print(filename_T)
                tumor_list.append(lst_T[0])
                img_t = mpimg.imread(os.path.join(root,filename_T))
                img_t.setflags(write=1)
                img_t[img_t == Counter(img_t.ravel()).most_common(1)[0][0]] = 0
                img_t[img_t > 0] = 2
                img_kt = img_t 
                KT_img.append(img_kt)
#                print(filename_K, 'does not exist')
                kidney_list.append([])
                #kidney_img.append(np.zeros((512,512),dtype=np.uint8))
        #        imgplot = plt.imshow(mpimg.imread(os.path.join(root,filename_T)))
        #        plt.show()            
            elif filename_K in files and not filename_T in files:
#                print(filename_T, 'does not exist')
                tumor_list.append([])
                #tumor_img.append(np.zeros((512,512),dtype=np.uint8))
#                print(filename_K)
                kidney_list.append(lst_K[0])
                img_k = mpimg.imread(os.path.join(root,filename_K))
                img_k.setflags(write=1)
                img_k[img_k == Counter(img_k.ravel()).most_common(1)[0][0]] = 0
                img_k[img_k > 0] = 1
                img_kt = img_k 
                KT_img.append(img_kt)
        #        imgplot = plt.imshow(mpimg.imread(os.path.join(root,filename_K)))
        #        plt.show()
            elif filename_K in files and filename_T in files:
#                print(filename_T)
                tumor_list.append(lst_T[0])
                img_t = mpimg.imread(os.path.join(root,filename_T))
                img_t.setflags(write=1)
                img_t[img_t == Counter(img_t.ravel()).most_common(1)[0][0]] = 0
                #tumor_img.append(img_t)
                img_t[img_t > 0] = 2
#                print(filename_K)
                kidney_list.append(lst_K[0])
                img_k = mpimg.imread(os.path.join(root,filename_K))
                img_k.setflags(write=1)
                img_k[img_k == Counter(img_k.ravel()).most_common(1)[0][0]] = 0
                #kidney_img.append(img_k) #-  mpimg.imread(os.path.join(root,filename_T)))
                img_k[img_k > 0] = 1
                
                img_kt = img_k + img_t
                img_kt[img_kt>2] = 1
                KT_img.append(img_kt)
        #        imgplot = plt.imshow(mpimg.imread(os.path.join(root,filename_T)))
        #        plt.show()     
        #        imgplot = plt.imshow(mpimg.imread(os.path.join(root,filename_K)))
        #        plt.show()
            else:     
                print('Segmented image not available') 
        
        CT_img_all.append(CT_img)
        subjname_CT = ''.join(['python_image_data/CT_img_subj',str(subj)])
        nsliceCT.append(np.shape(CT_img)[0])
        np.save(subjname_CT, CT_img)
        
        KT_img_all.append(KT_img)
        subjname_KT = ''.join(['python_image_data/KT_img_subj',str(subj)])
        nsliceKT.append(np.shape(CT_img)[0])
        np.save(subjname_KT, KT_img)

        CT_list_all.append(CT_list)
        tumor_list_all.append(tumor_list)
        kidney_list_all.append(kidney_list)
        
        subj += 1  

np.save('python_image_data/CT_list_all', CT_list_all)
np.save('python_image_data/tumor_list_all', tumor_list_all)
np.save('python_image_data/kidney_list_all', kidney_list_all)


        #######################################################################



#image = CT_img_all[0][2]        
#gray = cv2.GaussianBlur(image, (41, 41), 0)
#(minVal, maxVal, minLoc, maxLoc) = cv2.minMaxLoc(gray)
#shape = cv2.circle(image, maxLoc, 15, (0, 255, 0), 1)
#cv2.ellipse(image,maxLoc,(15,7),0,0,360,0,1)
#
##if shape.contains(maxLoc):
##    print('inside')
##else:
##    print('outside')    
#    
#
## display the results of our newly improved method
#cv2.imshow('rubost',image)
#cv2.waitKey(0)



        
        
        
          
