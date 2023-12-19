#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 11:19:55 2021

@author: matthewmartinez
Yi-Wei Chang Lab

Script to analyze general density plots. In this case (purpose of creation), analyzing density plots
through the central channel of the RSA
"""

import os
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
import re
import pandas as pd
import seaborn as sns

stream = os.popen('ls *.txt').read()
pointFiles = stream.split("\n"); pointFiles.pop()

allModels = {}

#Read AV data into a pandas dataframe, add all these to a dict (allModelsAV) with tomo ID as key
for x in pointFiles:
    with open(x,'r') as f:
        points = f.readlines()
    points = [re.split(r'\s+',line.strip()) for line in points]
    
    #Convert points into a Pandas dataframe
    df = pd.DataFrame({'Distance':[int(points[1][0])], 'Value':[float(points[1][1])]})
    
    count = 2
    for line in points[2:]:
        df = df.append(pd.DataFrame({'Distance':[int(points[count][0])], 'Value':[float(points[count][1])]}), ignore_index=True)
        count += 1
        
    name = os.path.basename(x) #replaces end of model file name with nothing, so nameTomo is the tomo ID
    allModels[name] = df
    
pixSize = 1.0616 #pixel size in nanometers

DensityProps = None
DensityProps = pd.DataFrame(columns = ['Name', 'RatioThird'])

pixelValuesList = []
for i,key in enumerate(allModels):
    pixelPoints = allModels[key]
    pixels = pixelPoints[['Distance','Value']].values
    
    Mems = argrelextrema(pixels[:,1], np.less)
    Mem1 = int(Mems[0][0])
    Mem2 = int(Mems[0][len(Mems[0])-1])
    
    PixBetween = [pixels[Mem1+1][1]]
    for k in range(Mem1+2,Mem2):
        PixBetween.append(pixels[k][1])
    
    #Create list of pixel values from 4 less than first membrane to 4 more than last membrane
    PixList = []
    if Mem1 >= 4:
        for k in range(Mem1-4,len(pixels)-1):
            PixList.append(pixels[k][1])
    else:
        for k in range(len(pixels)-1):
            PixList.append(pixels[k][1])
    
    if len(pixels)-1 - Mem2 >= 4:
        count = len(pixels)-1
        while count < Mem2+4:
            PixList.pop() 
            
    for i in range(len(PixList)):
        PixList[i] = 255-PixList[i]
    
    pixelValuesList.append(PixList)
            
    #Calculate average of center 1/3rd of plots between Central Channel
    Center = int((Mem1+Mem2)/2)
    PixThird = []
    for k in range(int(Center-(((Mem2-Mem1)/3)/2)),int(Center+(((Mem2-Mem1)/3)/2))):
        PixThird.append(pixels[k][1])
        
    RatioThird = np.mean(PixThird) / ((pixels[Mem1][1]+pixels[Mem2][1])/2)
    
    #Add all properties to dataframe
    DensityProps = DensityProps.append(pd.DataFrame({'Name':[key], 'RatioThird':[RatioThird]}), ignore_index=True)


keys_list = list(allModels)

Figure = plt.figure(i, figsize=(6,6))
ax = plt.axes([0.3,0.05,0.65,0.9])
ax.set_ylabel('Pixel Value',fontsize=60, labelpad=8)
ax.set_xlabel('Distance (nm)',fontsize=60, labelpad=8)
#ax1.set_axis_bgcolor(plt.cm.Blues(0.05))
ax.grid(False)
ax.tick_params('x', top=False, bottom=True, labelsize=16)
ax.tick_params('y', left=True, labelsize=16)
ax.set_ylim(0,255) #Find the max between AV and Rhop pixel values, set y lim to that +5
ax.set_xlim(0,100)
ax.plot(pixelValuesList[0])
ax.plot(pixelValuesList[1])
ax.plot(pixelValuesList[2])
