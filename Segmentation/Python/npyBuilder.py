# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 10:57:03 2018

@author: am43064
"""

import matplotlib.image as mpimg
from skimage import color
from skimage import io

img = []
for i in range(5):   
    filename = ''.join([str(i+1),'.jpg'])
    img_read = mpimg.imread(filename)
    img.append(color.rgb2gray(img_read))
    

np.save('testNPY', img)