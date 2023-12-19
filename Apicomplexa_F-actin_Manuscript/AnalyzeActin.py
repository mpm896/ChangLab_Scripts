#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 15:12:33 2022

@author: mattmartinez

Analyze cytosolic and pellicular F-actin (PCR-associated) in C. parvum

Object 1 = Cytosolic F-actin
Object 2 = Pellicular F-actin
Object 3 = Upper PCR to APR (Pt 1 = PCR, Pt 2 = APR)
    Ct 1 = Distance on one side of parasite
    Ct 2 = Distance on other side of parasite
    
ADDENDUM:
    If applicable, Object 4 = IMC collar density (APR) to actin filament distance for filaments channeled into the cytosol
    Each contour is a distance measurement. Pt 1 = actin filament, Pt 2 = IMC collar density/APR
"""

import numpy as np
import pandas as pd
import os
import sys
import glob
import re
import math
import scipy
import matplotlib.pyplot as plt
from matplotlib import cm
import random
import seaborn as sns
sns.set()

def calcDist(a,b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def calcTotalDist(a):
    cumDist=0
    for i in range(len(a)-1):
        dist = calcDist(a[i],a[i+1])
        cumDist += dist
    return cumDist

pointFiles = glob.glob("sm*/*_actin.txt", recursive=True) #These are for the "retracted conoid" tomograms
originalPointFiles = glob.glob("OriginalActinModels/*_actin.txt", recursive=True) #These are from the original actin segmentatio
allModels = {}
originalModels = {}

#Organize all model files into Pandas dataframe
for file in pointFiles:
    with open(file, 'r') as f:
        points = f.readlines()
    points = [re.split(r'\s+', line.strip()) for line in points]
    
    #convert points to Pandas dataframe
    df = pd.DataFrame({'Ob':[int(points[0][0])], 'Ct':[int(points[0][1])], 
    'X':[float(points[0][2])], 'Y':[float(points[0][3])], 'Z':[float(points[0][4])]})
    
    count = 1
    for line in points[1:]:
        df = df.append(pd.DataFrame({'Ob':[int(points[count][0])], 'Ct':[int(points[count][1])], 
        'X':[float(points[count][2])], 'Y':[float(points[count][3])], 
        'Z':[float(points[count][4])]}), ignore_index=True)        
        count += 1
    
    nameTomo = os.path.basename(file).replace('_actin.txt','')
    allModels[nameTomo]={'modelPoints':df}

#Original actin models
for file in originalPointFiles:
    with open(file, 'r') as f:
        points = f.readlines()
    points = [re.split(r'\s+', line.strip()) for line in points]
    
    #convert points to Pandas dataframe
    df = pd.DataFrame({'Ob':[int(points[0][0])], 'Ct':[int(points[0][1])], 
    'X':[float(points[0][2])], 'Y':[float(points[0][3])], 'Z':[float(points[0][4])]})
    
    count = 1
    for line in points[1:]:
        df = df.append(pd.DataFrame({'Ob':[int(points[count][0])], 'Ct':[int(points[count][1])], 
        'X':[float(points[count][2])], 'Y':[float(points[count][3])], 
        'Z':[float(points[count][4])]}), ignore_index=True)        
        count += 1
    
    nameTomo = os.path.basename(file).replace('_actin.txt','')
    originalModels[nameTomo]={'modelPoints':df}
    
ActinProps = None
ActinProps = pd.DataFrame(columns=['Name','NumCyto','NumPell','filRatio','LenCyto','LenPell','conDist1','conDist2','conDist3','conDist4'
                                'conDistAvg']) #LenCyto and LenPell will be lists

actinAPRdistances = []
for key in allModels:
    modelPoints = allModels[key]['modelPoints']
    
    cytoActin = modelPoints.loc[modelPoints['Ob'] == 1]
    pellActin = modelPoints.loc[modelPoints['Ob'] == 2]
    conoid = modelPoints.loc[modelPoints['Ob'] == 3] #PCR-APR measurement
    actinAPR = modelPoints.loc[modelPoints['Ob'] == 4] #Measurement from APR/IMC collar density to actin filament that is cytosolic
    
    if (len(conoid) != 8) or (conoid['Ct'].values[-1] != 4):
        print("***NOT ENOUGH CONOID POINTS***")
    numCyto = cytoActin['Ct'].values[-1]
    if len(pellActin['Ct'].values > 0):
        numPell = pellActin['Ct'].values[-1]
        filRatio = numCyto/(numCyto+numPell) #Ration of cytosolic to pellicular F-actin
        
        #Calculate length of each pellicular filament
        pellLengths = []
        for i in range(1,numPell+1):
            filPell = pellActin.loc[pellActin['Ct'] == i]
            filPellPoints = filPell[['X','Y','Z']].values
            lengthPell = calcTotalDist(filPellPoints)*1.06
            pellLengths.append(lengthPell)
    else:
        filRatio = 1
        numPell = 0
        pellLengths = []
            
    #Calculate length of each cytosolic filament
    cytoLengths = []
    for i in range(1,numCyto+1):
        filCyto = cytoActin.loc[cytoActin['Ct'] == i]
        filCytoPoints = filCyto[['X','Y','Z']].values
        lengthCyto = calcTotalDist(filCytoPoints)*1.06
        cytoLengths.append(lengthCyto)

    
   
        
    con1 = conoid.loc[conoid['Ct'] == 1]
    con1points = con1[['X','Y','Z']].values
    conDist1 = calcDist(con1points[0], con1points[1])*1.06
    
    con2 = conoid.loc[conoid['Ct'] == 2]
    con2points = con2[['X','Y','Z']].values
    conDist2 = calcDist(con2points[0], con2points[1])*1.06
    
    con3 = conoid.loc[conoid['Ct'] == 3]
    con3points = con3[['X','Y','Z']].values
    conDist3 = calcDist(con3points[0], con3points[1])*1.06
    
    con4 = conoid.loc[conoid['Ct'] == 4]
    con4points = con4[['X','Y','Z']].values
    conDist4 = calcDist(con4points[0], con4points[1])*1.06
    
    conLens = [conDist1,conDist2,conDist3,conDist4]
    conLens.sort() 
    conDistAvg = np.average(conLens)
    
    #Calculate the distance between the cytosolic actin filaments and the IMC collar
    if (len(actinAPR) > 0):
        for i in range(int(len(actinAPR)/2)):
            contour = actinAPR.loc[actinAPR['Ct'] == i+1]
            actinAPRpoints = contour[['X','Y','Z']].values
            actinAPRdist = calcDist(actinAPRpoints[0], actinAPRpoints[1])*1.06
            actinAPRdistances.append(actinAPRdist)
    else:
        actinAPRdist = float("NaN")
    
    ActinProps = ActinProps.append(pd.DataFrame({'Name':[key],'NumCyto':[numCyto],
    'NumPell':[numPell],'filRatio':[filRatio], 'LenCyto':[cytoLengths],'LenPell':[pellLengths],
    'conDist1':[conLens[0]], 'conDist2':[conLens[1]],'conDist3':[conLens[2]], 'conDist4':[conLens[3]], 'conDistAvg':[conDistAvg], 'actinAPR':[actinAPRdist]}),ignore_index=True)

cytoActin = None
pellActin = None
conoid = None
actinAPR = None

for key in originalModels:
    modelPoints = originalModels[key]['modelPoints']
    
    pellActin = modelPoints.loc[modelPoints['Ob'] == 1]
    cytoActin = modelPoints.loc[modelPoints['Ob'] == 2]
    conoid = modelPoints.loc[modelPoints['Ob'] == 3] #PCR-APR measurement
    actinAPR = modelPoints.loc[modelPoints['Ob'] == 4] #Measurement from APR/IMC collar density to actin filament that is cytosolic
    
    numPell = pellActin['Ct'].values[-1]
    if len(cytoActin['Ct'].values > 0):
        numCyto = cytoActin['Ct'].values[-1]
        filRatio = numCyto/(numCyto+numPell) #Ration of cytosolic to pellicular F-actin
        
        #Calculate length of each pellicular filament
        cytoLengths = []
        for i in range(1,numCyto+1):
            filCyto = cytoActin.loc[cytoActin['Ct'] == i]
            filCytoPoints = filCyto[['X','Y','Z']].values
            lengthCyto = calcTotalDist(filCytoPoints)*1.06
            cytoLengths.append(lengthCyto)
    else:
        filRatio = 0
        numCyto = 0
        cytoLengths = []
            
    #Calculate length of each cytosolic filament
    pellLengths = []
    for i in range(1,numPell+1):
        filPell = pellActin.loc[pellActin['Ct'] == i]
        filPellPoints = filPell[['X','Y','Z']].values
        lengthPell = calcTotalDist(filPellPoints)*1.06
        #if (key == 'sm2019-09-07-4'):
            #lengthPell /= 10
        pellLengths.append(lengthPell)
        
    con1 = conoid.loc[conoid['Ct'] == 1]
    con1points = con1[['X','Y','Z']].values
    conDist1 = calcDist(con1points[0], con1points[1])*1.06
    
    con2 = conoid.loc[conoid['Ct'] == 2]
    con2points = con2[['X','Y','Z']].values
    conDist2 = calcDist(con2points[0], con2points[1])*1.06
    
    con3 = conoid.loc[conoid['Ct'] == 3]
    con3points = con3[['X','Y','Z']].values
    conDist3 = calcDist(con3points[0], con3points[1])*1.06
    
    con4 = conoid.loc[conoid['Ct'] == 4]
    con4points = con4[['X','Y','Z']].values
    conDist4 = calcDist(con4points[0], con4points[1])*1.06
    
    conLens = [conDist1,conDist2,conDist3,conDist4]
    conLens.sort() 
    conDistAvg = np.average(conLens)
    
    #Calculate the distance between the cytosolic actin filaments and the IMC collar
    if (len(actinAPR) > 0):
        for i in range(int(len(actinAPR)/2)):
            contour = actinAPR.loc[actinAPR['Ct'] == i+1]
            actinAPRpoints = contour[['X','Y','Z']].values
            actinAPRdist = calcDist(actinAPRpoints[0], actinAPRpoints[1])*1.06
            actinAPRdistances.append(actinAPRdist)
    else:
        actinAPRdist = float("NaN")
    
    ActinProps = ActinProps.append(pd.DataFrame({'Name':[key],'NumCyto':[numCyto],
    'NumPell':[numPell],'filRatio':[filRatio], 'LenCyto':[cytoLengths],'LenPell':[pellLengths],
    'conDist1':[conLens[0]], 'conDist2':[conLens[1]],'conDist3':[conLens[2]], 'conDist4':[conLens[3]], 'conDistAvg':[conDistAvg], 'actinAPR':[actinAPRdist]}),ignore_index=True)

#Calculate total number of filaments in each class, and fraction of cytosolic filaments in which we could measure the distance to the IMC collar
totalCyto = ActinProps['NumCyto'].sum()
totalPell = ActinProps['NumPell'].sum()
fracCyto = len(actinAPRdistances)/totalCyto

#Plot of cytosolic F-actin lengths    
Figure1 = plt.figure(1, figsize=(4,4))
ax1 = plt.axes([0.3,0.05,0.65,0.9])
ax1.set_ylabel('PCR-Actin lengths (nm)',fontsize=30, labelpad=8)
ax1.set_facecolor(cm.Blues(0.05))
ax1.grid(False)
ax1.tick_params('x', top=False, bottom=False, labelsize=0)
ax1.tick_params('y', right=False, labelsize=16)

cytoActinLengths = np.concatenate(ActinProps['LenCyto'],axis=0)
pellicularActinLengths = np.concatenate(ActinProps['LenPell'],axis=0)

ax1.set_ylim(-25,max(cytoActinLengths)+15)
ax1.set_xlim(-0.3,1.3)
ax1_1 = sns.violinplot(data=[cytoActinLengths,pellicularActinLengths],color=cm.Blues(0.2))
ax1_2 = sns.swarmplot(data=[cytoActinLengths,pellicularActinLengths],s=4.5,color='b',alpha=0.9)
#plt.savefig('/ChangLab1-hd2/matt/Cryptosporidium/ApicalRing/RetractedConoid/ActinLengths.png', dpi=600, bbox_inches='tight')

#Plot ratio of filaments vs conoid extrusion length
Figure2 = plt.figure(2, figsize=(4,4))
ax2 = plt.axes([0.3,0.05,0.65,0.9])
ax2.set_ylabel('Proportion of F-actin in cytosol',fontsize=20,labelpad=8)
ax2.set_xlabel('Distance from upper PCR to APR (nm)',fontsize=15,labelpad=8)
ax2.set_facecolor(cm.Blues(0.05))
ax2.grid(False)
filRatios = ActinProps['filRatio'].values

ax2.plot(ActinProps['conDistAvg'], filRatios, 'o', mfc='b', mec='k', markersize=8, 
         alpha=0.9)

x = ActinProps['conDistAvg'].values
#theta = np.polyfit(x, filRatios, 1)
#y_line = theta[1] + theta[0] * x
#plt.plot(x, y_line, 'r')
ax2.plot(np.unique(x), np.poly1d(np.polyfit(x, filRatios.tolist(), 1))(np.unique(x)),color='r',linestyle='--')
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, filRatios.tolist())
rsquared = r_value**2
ax2.set_ylim(-0.1,1.1)

#plt.savefig('/ChangLab1-hd2/matt/Cryptosporidium/ApicalRing/RetractedConoid/ActinLocation.png', dpi=600, bbox_inches='tight')

Figure3 = plt.figure(3, figsize=(4,4))
ax3 = plt.axes([0.3,0.05,0.65,0.9])

ax3.plot(ActinProps['conDist1'], filRatios, 'o', mfc='b', mec='k', markersize=8, 
         alpha=0.4)

Figure4 = plt.figure(4, figsize=(4,4))
ax4 = plt.axes([0.3,0.05,0.65,0.9])

ax4.plot(ActinProps['conDist4'], filRatios, 'o', mfc='b', mec='k', markersize=8, 
         alpha=0.4)

#Plot IMC collar - actin filament distances
Figure5 = plt.figure(5, figsize=(4,4))
ax5 = plt.axes([0.3,0.05,0.65,0.9])
ax5.set_ylabel('APR - F-actin distance (nm)',fontsize=30, labelpad=8)
#ax5.set_axis_bgcolor(plt.cm.Blues(0.05))
ax5.grid(False)
ax5.tick_params('x', top=False, bottom=False, labelsize=0)
ax5.tick_params('y', right=False, labelsize=16)
ax5.set_ylim(-1,max(actinAPRdistances)+5)
ax1.set_xlim(-0.3,0.3)
meanActinAPRdist = np.mean(actinAPRdistances)
stdActinAPRdist = np.std(actinAPRdistances)
#ax5.plot([-0.53,0.53],[meanActinAPRdist,meanActinAPRdist],color='k',ls='-',lw=5,alpha=0.7)
ax6 = sns.violinplot(data=actinAPRdistances, color='b')
plt.setp(ax6.collections, alpha=0.3)
ax5 = sns.swarmplot(data=actinAPRdistances,s=12,color='b')
