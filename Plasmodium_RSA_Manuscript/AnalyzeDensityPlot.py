#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 17:14:45 2021

@author: Matt Martinez

This script will analysze a density plot profile across a straight line. Will analyze:
    1. Stats on the density between membranes
    2. Widths of membranes. 
    3. Unclear yet all that it will analyze. Stay tuned. 
    
INPUT: 2 column text file. Column 1 is pixel number, column 2 is pixel value

To get the proper files for input, my workflow:
    1. Take image of slice of tomogram, save as TIFF
    2. Load into Image J (or FIJI). Convert to 32-bit grayscale: -->This is unnecessary
        Image -> type -> 32-bit
    3. Apply Gaussian filter:
        Process -> filter -> Gaussian (2 pixels is fine)
    4. Make a line across the density you want to analyze. In my case, drew line across the AV and rhoptry
        Set line width to something greater than 1. This will average x pixels orthoganally to that point in the line (reduces noise)
    5. After making line and changing width, go to analyze -> plot profile
    6. Set to live and adjust as necessary. Output to text file:
        List -> file, save OR
        Data -> save (as txt)
        
Adjust pixel size in line 82 to that of your tomograms to be analyzed
"""

import os
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
import re
import pandas as pd
import seaborn as sns

sns.set()

streamAV = os.popen('ls *_AV.txt').read()
pointFiles_AV = streamAV.split("\n"); pointFiles_AV.pop()
streamRhop = os.popen('ls *_rhoptry.txt').read()
pointFiles_Rhop = streamRhop.split("\n"); pointFiles_Rhop.pop()

allModelsAV = {}
allModelsRhop = {}

#Read AV data into a pandas dataframe, add all these to a dict (allModelsAV) with tomo ID as key
for x in pointFiles_AV:
    with open(x,'r') as f:
        points = f.readlines()
    points = [re.split(r'\s+',line.strip()) for line in points]
    
    #Convert points into a Pandas dataframe
    df = pd.DataFrame({'Distance':[int(points[1][0])], 'Value':[float(points[1][1])]})
    
    count = 2
    for line in points[2:]:
        df = df.append(pd.DataFrame({'Distance':[int(points[count][0])], 'Value':[float(points[count][1])]}), ignore_index=True)
        count += 1
        
    nameTomo = os.path.basename(x).replace('_AV.txt','') #replaces end of model file name with nothing, so nameTomo is the tomo ID
    allModelsAV[nameTomo] = df
    
#Read rhoptry data into a pandas dataframe, add all these to a dict (allModelsRhop) with tomo ID as key
for x in pointFiles_Rhop:
    with open(x,'r') as f:
        points = f.readlines()
    points = [re.split(r'\s+',line.strip()) for line in points]
    
    #Convert points into a Pandas dataframe
    df = pd.DataFrame({'Distance':[int(points[1][0])], 'Value':[float(points[1][1])]})
    
    count = 2
    for line in points[2:]:
        df = df.append(pd.DataFrame({'Distance':[int(points[count][0])], 'Value':[float(points[count][1])]}), ignore_index=True)
        count += 1
        
    nameTomo = os.path.basename(x).replace('_rhoptry.txt','') #replaces end of model file name with nothing, so nameTomo is the tomo ID
    allModelsRhop[nameTomo] = df
    
pixSize = 1.0616 #pixel size in nanometers

#Recomplile density plot properties as follows
DensityProps = None
DensityProps = pd.DataFrame(columns = ['Name', 'AvRatio1', 'AvRatio2', 'AvRatioAvg','AvRatioThird', 'AvDiff1', 'AvDiff2', 'AvDiffAvg', 'RhopRatio1', 'RhopRatio2', 'RhopRatioAvg','RhopRatioThird', 'RhopDiff1', 'RhopDiff2', 'RhopDiffAvg'])

for i,key in enumerate(allModelsAV):
    pixelPointsAV = allModelsAV[key]
    pixelPointsRhoptry = allModelsRhop[key]
    
    pixelsAV = pixelPointsAV[['Distance','Value']].values
    pixelsRhop = pixelPointsRhoptry[['Distance','Value']].values
    
    
    
    #Invert pixelsAV and pixelsRhop so that 0 is less dense and 255 is more dense
    for i in range(len(pixelsAV)):
        pixelsAV[i][1] = 255-pixelsAV[i][1]
        
    for i in range(len(pixelsRhop)):
        pixelsRhop[i][1] = 255-pixelsRhop[i][1]
    
    #Fine first and last maxima that represent AV membrane points
    AvMems = argrelextrema(pixelsAV[:,1], np.greater)
    AvMem1 = int(AvMems[0][0])
    AvMem2 = int(AvMems[0][len(AvMems[0])-1])
    
    AvBetween = [pixelsAV[AvMem1+1][1]]
    for k in range(AvMem1+2,AvMem2):
        AvBetween.append(pixelsAV[k][1])
    
    #Create list of AV pixel values from 4 less than first membrane to 4 more than last membrane
    AvList = []
    if AvMem1 >= 4:
        for k in range(AvMem1-4,len(pixelsAV)-1):
            AvList.append(pixelsAV[k][1])
    else:
        for k in range(len(pixelsAV)-1):
            AvList.append(pixelsAV[k][1])
    
    if len(pixelsAV)-1 - AvMem2 >= 4:
        count = len(pixelsAV)-1
        while count < AvMem2+4:
            AvList.pop()        
    
    #Fine first and last minima that represent AV membrane points
    RhopMems = argrelextrema(pixelsRhop[:,1], np.greater)
    RhopMem1 = int(RhopMems[0][0])
    RhopMem2 = int(RhopMems[0][len(RhopMems[0])-1])
    
    RhopBetween = [pixelsRhop[RhopMem1+1][1]]
    for k in range(RhopMem1+2,RhopMem2):
        RhopBetween.append(pixelsRhop[k][1])
        
    #Create list of AV pixel values from 4 less than the first membrane to 4 more than the last membrane
    RhopList = []
    if RhopMem1 >= 4:
        for k in range(RhopMem1-4,len(pixelsRhop)-1):
            RhopList.append(pixelsRhop[k][1])
    else:
        for k in range(len(pixelsRhop)-1):
            RhopList.append(pixelsRhop[k][1])
    
    if len(pixelsRhop)-1 - RhopMem2 >= 4:
        count = len(pixelsRhop)-1
        while count < RhopMem2+4:
            RhopList.pop()    

    #Calculate average of center 1/3rd of Av/Rhoptry plots between membranes
    AvCenter = int((AvMem1+AvMem2)/2)
    AvThird = []
    for k in range(int(AvCenter-(((AvMem2-AvMem1)/3)/2)),int(AvCenter+(((AvMem2-AvMem1)/3)/2))):
        AvThird.append(pixelsAV[k][1])

    RhopCenter = int((RhopMem1+RhopMem2)/2)  
    RhopThird = []
    for k in range(int(RhopCenter-(((RhopMem2-RhopMem1)/3)/2)),int(RhopCenter+(((RhopMem2-RhopMem1)/3)/2))):
        RhopThird.append(pixelsRhop[k][1])
    
    
    #Calculate ratios from each membrane to maximum and from membrane average to maximum
    AvMax = max(AvBetween)
    AvRatio1 = AvMax/pixelsAV[AvMem1][1]
    AvDiff1 = AvMax-pixelsAV[AvMem1][1]
    AvRatio2 = AvMax/pixelsAV[AvMem2][1]
    AvDiff2 = AvMax-pixelsAV[AvMem2][1]
    AvRatioAvg = AvMax / ((pixelsAV[AvMem1][1]+pixelsAV[AvMem2][1])/2)
    AvDiffAvg = AvMax - ((pixelsAV[AvMem1][1]+pixelsAV[AvMem2][1])/2)
    AvRatioThird = np.mean(AvThird) / ((pixelsAV[AvMem1][1]+pixelsAV[AvMem2][1])/2)
    
    RhopMax = max(RhopBetween)
    RhopRatio1 = RhopMax/pixelsRhop[RhopMem1][1]
    RhopDiff1 = RhopMax-pixelsRhop[RhopMem1][1]
    RhopRatio2 = RhopMax/pixelsRhop[RhopMem2][1]
    RhopDiff2 = RhopMax-pixelsRhop[RhopMem2][1]
    RhopRatioAvg = RhopMax / ((pixelsRhop[RhopMem1][1]+pixelsRhop[RhopMem2][1])/2)
    RhopDiffAvg = RhopMax - ((pixelsRhop[RhopMem1][1]+pixelsRhop[RhopMem2][1])/2)
    RhopRatioThird = np.mean(RhopThird) / ((pixelsRhop[RhopMem1][1]+pixelsRhop[RhopMem2][1])/2)
    

    
    Figure = plt.figure(i, figsize=(4,4))
    ax = plt.axes([0.3,0.05,0.65,0.9])
    ax.set_ylabel('Pixel Value',fontsize=30, labelpad=8)
    ax.set_xlabel('Distance (nm)',fontsize=30, labelpad=8)
    #ax1.set_axis_bgcolor(plt.cm.Blues(0.05))
    ax.grid(False)
    ax.tick_params('x', top=False, bottom=True, labelsize=16)
    ax.tick_params('y', left=True, labelsize=16)
    ax.set_ylim(min(min(AvList),min(RhopList))-20,max(max(AvList),max(RhopList))+20) #Find the max between AV and Rhop pixel values, set y lim to that +5
    ax.set_xlim(0,len(RhopList))
    ax.plot(range(len(AvList)), AvList)
    ax.plot(range(len(RhopList)), RhopList)
    #plt.gca().invert_yaxis()
    
    #plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/densityPlot_'+key+'.png', bbox_inches="tight") #bbox_inches made it so the saved figure wasn't cut off
    
    #Print text: Av and rhoptry membranes, and average density between membranes
    print("Average Av density in center 1/3rd is ",np.mean(AvThird))
    print("Average rhoptry density in center 1/3rd is ",np.mean(RhopThird))
    print("Av ratio with average membrane: ",AvRatioAvg)
    print("Rhop ratio with average membrane: ",RhopRatioAvg)
    
    #Add all properties to dataframe
    DensityProps = DensityProps.append(pd.DataFrame({'Name':[key], 'AvRatio1':[AvRatio1], 'AvRatio2':[AvRatio2], 'AvRatioAvg':[AvRatioAvg], 'AvRatioThird':[AvRatioThird], 'AvDiff1':[AvDiff1], 'AvDiff2':[AvDiff2], 'AvDiffAvg':[AvDiffAvg], 'RhopRatio1':[RhopRatio1], 'RhopRatio2':[RhopRatio2], 'RhopRatioAvg':[RhopRatioAvg], 'RhopRatioThird':[RhopRatioThird], 'RhopDiff1':[RhopDiff1], 'RhopDiff2':[RhopDiff2], 'RhopDiffAvg':[RhopDiffAvg]}), ignore_index=True)

    #For source data
    print("\n",key)
    print("\nAV Ratios: ", DensityProps['AvRatioThird'])
    print("\nRhoptry Ratios: ", DensityProps['RhopRatioThird'])
    #print("\nAV Pixel Values: ", AvList)
    #print("\nRhoptry Pixel Values: ", RhopList)
'''       
Figure1 = plt.figure(15, figsize=(6,6))
ax15 = plt.axes([0.3,0.05,0.65,0.9])
ax15.set_ylabel('Ratio',fontsize=30, labelpad=8)
plt.title("Ratio of peak intensity to membrane, AV vs. Rhoptry",fontsize=16)
#ax1.set_axis_bgcolor(plt.cm.Blues(0.05))
ax15.grid(False)
ax15.tick_params('x', top=False, bottom=True, labelsize=16)
ax15.tick_params('y', left=True, labelsize=16)
meanAVratio = np.mean(DensityProps['AvRatioAvg'].values)
meanRhopratio = np.mean(DensityProps['RhopRatioAvg'].values)
#ax15.set_xlim(-10,10)
#ax15.plot([-0.53,0.53],[meanAVratio,meanAVratio], color='b',ls='--',lw=5,alpha=0.7)
#ax15.plot([0.5,1],[meanRhopratio,meanRhopratio],color='darkorange',ls='--',lw=5,alpha=0.7)
ax15 = sns.swarmplot(data=[DensityProps['AvRatioAvg'],DensityProps['RhopRatioAvg']],s=8,alpha=0.9)
sns.boxplot(data=[DensityProps['AvRatioAvg'],DensityProps['RhopRatioAvg']], color='w', showmeans=True, meanline=True, meanprops={"lw":2,"color":"black"})
#markers, caps, bars = ax1.errorbar(x=0,y=meanAvDist,yerr=stdAvDist,elinewidth=5,ecolor='r',capsize=30,capthick=5)  
#ax15.plot(DensityProps['AvRatio1'])
#bar.set_alpha(0.99) for bar in bars]
#[cap.set_alpha(0.99) for cap in caps] 
plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/DensityRatio.png')  


Figure2 = plt.figure(16,figsize=(6,6))
ax16 = plt.axes([0.3,0.05,0.65,0.9])
ax16.set_ylabel('Difference (pixel values)',fontsize=30, labelpad=8)
plt.title("Difference between peak intensity and membrane, AV vs. Rhoptry", fontsize=16)
ax16.grid(False)
ax16.tick_params('x', top=False, bottom=True, labelsize=16)
ax16.tick_params('y', left=True, labelsize=16)
meanAvDiff = np.mean(DensityProps['AvDiffAvg'].values)
meanRhopDiff = np.mean(DensityProps['RhopDiffAvg'].values)
ax16 = sns.swarmplot(data=[DensityProps['AvDiffAvg'],DensityProps['RhopDiffAvg']],s=8,alpha=0.9)
sns.boxplot(data=[DensityProps['AvDiffAvg'],DensityProps['RhopDiffAvg']], color='w', showmeans=True, meanline=True, meanprops={"lw":2,"color":"black"})

plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/DensityDifference.png')  
'''

Figure4 = plt.figure(17,figsize=(4,6))
ax17 = plt.axes([0.3,0.05,0.8,0.6])
ax17.set_ylabel('Ratio',fontsize=30, labelpad=8)
#plt.title("Ratio of internal density to membrane density, AV vs rhoptry", fontsize=16)
ax17.grid(False)
ax17.tick_params('x', top=False, bottom=False, labelsize=0)
ax17.tick_params('y', left=True, labelsize=16)
clrs = ['#5975a4e5','#cc8863e6']
ax17 = sns.swarmplot(data=[DensityProps['AvRatioThird'],DensityProps['RhopRatioThird']],s=12,palette=clrs,alpha=0.9)
sns.boxplot(data=[DensityProps['AvRatioThird'],DensityProps['RhopRatioThird']], color='w', showmeans=True, meanline=True, meanprops={"lw":2,"color":"black"})

#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/AverageThirdDensityRatio.png')  

    
    

    