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

streamAV = os.popen('ls AV/Analysis/*_AV.txt').read()
pointFiles_AV = streamAV.split("\n"); pointFiles_AV.pop()

streamPre = os.popen('ls Prefusion/Analysis/*_Prefusion.txt').read()
pointFiles_Pre = streamPre.split("\n"); pointFiles_Pre.pop()

streamFus = os.popen('ls Fusion/Analysis/*_Fusion.txt').read()
pointFiles_Fus = streamFus.split("\n"); pointFiles_Fus.pop()

allModelsAV = {}
allModelsPre = {}
allModelsFus = {}

#Read AV rhoptry data into a pandas dataframe, add all these to a dict (allModelsAV) with tomo ID as key
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
    
#Read prefusion rhoptry data into a pandas dataframe, add all these to a dict (allModelsRhop) with tomo ID as key
for x in pointFiles_Pre:
    with open(x,'r') as f:
        points = f.readlines()
    points = [re.split(r'\s+',line.strip()) for line in points]
    
    #Convert points into a Pandas dataframe
    df = pd.DataFrame({'Distance':[int(points[1][0])], 'Value':[float(points[1][1])]})
    
    count = 2
    for line in points[2:]:
        df = df.append(pd.DataFrame({'Distance':[int(points[count][0])], 'Value':[float(points[count][1])]}), ignore_index=True)
        count += 1
        
    nameTomo = os.path.basename(x).replace('_Prefusion.txt','') #replaces end of model file name with nothing, so nameTomo is the tomo ID
    allModelsPre[nameTomo] = df
    
#Read prefusion rhoptry data into a pandas dataframe, add all these to a dict (allModelsRhop) with tomo ID as key
for x in pointFiles_Fus:
    with open(x,'r') as f:
        points = f.readlines()
    points = [re.split(r'\s+',line.strip()) for line in points]
    
    #Convert points into a Pandas dataframe
    df = pd.DataFrame({'Distance':[int(points[1][0])], 'Value':[float(points[1][1])]})
    
    count = 2
    for line in points[2:]:
        df = df.append(pd.DataFrame({'Distance':[int(points[count][0])], 'Value':[float(points[count][1])]}), ignore_index=True)
        count += 1
        
    nameTomo = os.path.basename(x).replace('_Fusion.txt','') #replaces end of model file name with nothing, so nameTomo is the tomo ID
    allModelsFus[nameTomo] = df
    
pixSize = 1.0616 #pixel size in nanometers

#Recomplile AV rhoptry density plot properties as follows
AvDensityProps = None
AvDensityProps = pd.DataFrame(columns = ['Name', 'Mem1', 'Mem2', 'MemAvg', 'RhopAvg', 'Ratio'])

#Recomplile Prefusion rhoptry density plot properties as follows
PreDensityProps = None
PreDensityProps = pd.DataFrame(columns = ['Name', 'Mem1', 'Mem2', 'MemAvg', 'RhopAvg', 'Ratio'])

#Recomplile Fusion rhoptry density plot properties as follows
FusDensityProps = None
FusDensityProps = pd.DataFrame(columns = ['Name', 'Mem1', 'Mem2', 'MemAvg', 'RhopAvg', 'Ratio'])

for i,key in enumerate(allModelsAV):
    pixelPointsAV = allModelsAV[key]
    pixelsAV = pixelPointsAV[['Distance','Value']].values
    
    #Invert pixelsAV and pixelsRhop so that 0 is less dense and 255 is more dense
    for i in range(len(pixelsAV)):
        pixelsAV[i][1] = 255-pixelsAV[i][1]
    
    #Fine first and last minima that represent AV membrane points
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
    
    
    #Calculate average of center 1/3rd of Av/Rhoptry plots between membranes        
    AvCenter = int((AvMem1+AvMem2)/2)
    AvThird = []
    for k in range(int(AvCenter-(((AvMem2-AvMem1)/3)/2)),int(AvCenter+(((AvMem2-AvMem1)/3)/2))):
        AvThird.append(pixelsAV[k][1])

    AvRatio = np.mean(AvThird) / ((pixelsAV[AvMem1][1]+pixelsAV[AvMem2][1])/2)
    MemAvg = (pixelsAV[AvMem1][1]+pixelsAV[AvMem2][1])/2
    
    Figure = plt.figure(i, figsize=(6,6))
    ax = plt.axes([0.3,0.05,0.65,0.9])
    ax.set_ylabel('Pixel Value',fontsize=20, labelpad=8)
    ax.set_xlabel('Distance (nm)',fontsize=20, labelpad=8)
    #ax1.set_axis_bgcolor(plt.cm.Blues(0.05))
    ax.grid(False)
    ax.tick_params('x', top=False, bottom=True, labelsize=16)
    ax.tick_params('y', left=True, labelsize=16)
    ax.set_ylim(0,255) #Find the max between AV and Rhop pixel values, set y lim to that +5
    ax.set_xlim(0,len(pixelsAV)+2)
    ax.plot(range(len(AvList)), AvList)
    
    #Add all properties to dataframe
    AvDensityProps = AvDensityProps.append(pd.DataFrame({'Name':[key], 'Mem1':[pixelsAV[AvMem1][1]], 'Mem2':[pixelsAV[AvMem2][1]], 'MemAvg':[MemAvg], 'RhopAvg':[np.mean(AvThird)], 'Ratio':[AvRatio]}), ignore_index=True)


for i,key in enumerate(allModelsPre):
    pixelPointsPre = allModelsPre[key]
    pixelsPre = pixelPointsPre[['Distance','Value']].values
    
    #Invert pixelsAV and pixelsRhop so that 0 is less dense and 255 is more dense
    for i in range(len(pixelsPre)):
        pixelsPre[i][1] = 255-pixelsPre[i][1]
    
    #Fine first and last minima that represent AV membrane points
    Mems = argrelextrema(pixelsPre[:,1], np.greater)
    Mem1 = int(Mems[0][0])
    Mem2 = int(Mems[0][len(Mems[0])-1])
    
    MemBetween = [pixelsPre[Mem1+1][1]]
    for k in range(Mem1+2,Mem2):
        MemBetween.append(pixelsPre[k][1])
    
    #Create list of AV pixel values from 4 less than first membrane to 4 more than last membrane
    PreList = []
    if Mem1 >= 4:
        for k in range(Mem1-4,len(pixelsPre)-1):
            PreList.append(pixelsPre[k][1])
    else:
        for k in range(len(pixelsPre)-1):
            PreList.append(pixelsPre[k][1])
    
    if len(pixelsPre)-1 - Mem2 >= 4:
        count = len(pixelsPre)-1
        while count < Mem2+4:
            PreList.pop()   
    
    for i in range(len(PreList)):
        PreList[i] = 255-PreList[i]
    
    #Calculate average of center 1/3rd of Av/Rhoptry plots between membranes        
    PreCenter = int((Mem1+Mem2)/2)
    PreThird = []
    for k in range(int(PreCenter-(((Mem2-Mem1)/3)/2)),int(PreCenter+(((Mem2-Mem1)/3)/2))):
        PreThird.append(pixelsPre[k][1])

    Ratio = np.mean(PreThird) / ((pixelsPre[Mem1][1]+pixelsPre[Mem2][1])/2)
    MemAvg = (pixelsPre[Mem1][1]+pixelsPre[Mem2][1])/2
    
    Figure = plt.figure(i+10, figsize=(6,6))
    ax = plt.axes([0.3,0.05,0.65,0.9])
    ax.set_ylabel('Pixel Value',fontsize=20, labelpad=8)
    ax.set_xlabel('Distance (nm)',fontsize=20, labelpad=8)
    #ax1.set_axis_bgcolor(plt.cm.Blues(0.05))
    ax.grid(False)
    ax.tick_params('x', top=False, bottom=True, labelsize=16)
    ax.tick_params('y', left=True, labelsize=16)
    ax.set_ylim(0,255) #Find the max between AV and Rhop pixel values, set y lim to that +5
    ax.set_xlim(0,len(pixelsPre)+2)
    ax.plot(range(len(PreList)), PreList)
    
    #Add all properties to dataframe
    PreDensityProps = PreDensityProps.append(pd.DataFrame({'Name':[key], 'Mem1':[pixelsPre[Mem1][1]], 'Mem2':[pixelsPre[Mem2][1]], 'MemAvg':[MemAvg], 'RhopAvg':[np.mean(PreThird)], 'Ratio':[Ratio]}), ignore_index=True)


for i,key in enumerate(allModelsFus):
    pixelPointsFus = allModelsFus[key]
    pixelsFus = pixelPointsFus[['Distance','Value']].values
    
    #Invert pixelsAV and pixelsRhop so that 0 is less dense and 255 is more dense
    for i in range(len(pixelsFus)):
        pixelsFus[i][1] = 255-pixelsFus[i][1]
        
    #Fine first and last minima that represent AV membrane points
    Mems = argrelextrema(pixelsFus[:,1], np.greater)
    Mem1 = int(Mems[0][0])
    Mem2 = int(Mems[0][len(Mems[0])-1])
    
    MemBetween = [pixelsFus[Mem1+1][1]]
    for k in range(Mem1+2,Mem2):
        MemBetween.append(pixelsFus[k][1])
    
    #Create list of AV pixel values from 4 less than first membrane to 4 more than last membrane
    FusList = []
    if Mem1 >= 4:
        for k in range(Mem1-4,len(pixelsFus)-1):
            FusList.append(pixelsFus[k][1])
    else:
        for k in range(len(pixelsFus)-1):
            FusList.append(pixelsFus[k][1])
    
    if len(pixelsFus)-1 - Mem2 >= 4:
        count = len(pixelsFus)-1
        while count < Mem2+4:
            FusList.pop()
            
    for i in range(len(FusList)):
        FusList[i] = 255-FusList[i]
    
    #Calculate average of center 1/3rd of Av/Rhoptry plots between membranes        
    FusCenter = int((Mem1+Mem2)/2)
    FusThird = []
    for k in range(int(FusCenter-(((Mem2-Mem1)/3)/2)),int(FusCenter+(((Mem2-Mem1)/3)/2))):
        FusThird.append(pixelsFus[k][1])

    Ratio = np.mean(FusThird) / ((pixelsFus[Mem1][1]+pixelsFus[Mem2][1])/2)
    MemAvg = (pixelsFus[Mem1][1]+pixelsFus[Mem2][1])/2
    
    Figure = plt.figure(i+19, figsize=(6,6))
    ax = plt.axes([0.3,0.05,0.65,0.9])
    ax.set_ylabel('Pixel Value',fontsize=20, labelpad=8)
    ax.set_xlabel('Distance (nm)',fontsize=20, labelpad=8)
    #ax1.set_axis_bgcolor(plt.cm.Blues(0.05))
    ax.grid(False)
    ax.tick_params('x', top=False, bottom=True, labelsize=16)
    ax.tick_params('y', left=True, labelsize=16)
    ax.set_ylim(0,255) #Find the max between AV and Rhop pixel values, set y lim to that +5
    ax.set_xlim(0,len(pixelsFus)+2)
    ax.plot(range(len(FusList)), FusList)
    
    #Add all properties to dataframe
    FusDensityProps = FusDensityProps.append(pd.DataFrame({'Name':[key], 'Mem1':[pixelsFus[Mem1][1]], 'Mem2':[pixelsFus[Mem2][1]], 'MemAvg':[MemAvg], 'RhopAvg':[np.mean(FusThird)], 'Ratio':[Ratio]}), ignore_index=True)


Figure29 = plt.figure(29,figsize=(4,4))
ax17 = plt.axes([0.3,0.05,0.8,0.6])
ax17.set_ylabel('Ratio',fontsize=30, labelpad=8)
#plt.title("Ratio of internal density to membrane density, AV vs rhoptry", fontsize=16)
ax17.grid(False)
ax17.tick_params('x', top=False, bottom=False, labelsize=0)
#locs, labels = plt.xticks()
#plt.xticks([0.4, 0.8], ['AV', 'Prefusion', 'Fusion'])
ax17.set_xlim(-0.3,0.3)
ax17.tick_params('y', left=True, labelsize=16)
ax17 = sns.swarmplot(data=[AvDensityProps['Ratio'],PreDensityProps['Ratio'],FusDensityProps['Ratio']],s=12,color='b',alpha=0.9)
sns.boxplot(data=[AvDensityProps['Ratio'],PreDensityProps['Ratio'],FusDensityProps['Ratio']], color='w', showmeans=True, meanline=True, meanprops={"lw":2,"color":"black"})

