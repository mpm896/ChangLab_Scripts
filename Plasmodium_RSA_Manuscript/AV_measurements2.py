#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
9/16/2021
Matthew Martinez
Chang Lab

Extract values from apicomplexan AV measurements.
OBJECT 1, CONTOUR 1: 2 points, at PM apex and closest point on the AV
OBJECT 2, CONTOUR 1: 2 points along long axis of the AV
OBJECT 2, CONTOUR 2: 2 points along the short axis of the AV
OBJECT 3, CONTOUR 1: 1 point at the tip of the rhoptry closest to the long axis of the AV
OBJECT 3, CONTOUR 2: 2 points, one at the tip of the other docked rhoptry and 1 at the closest point on the AV

All objects are open contour type.
This script will measure these values:
    1. AV long and short axis lengths and eccentricity
    2. Distance from PM apex to closest point of the AV
    3. Angle at which the AV is docked at the RSA
    4. Distance of both rhoptries (if docked) to the AV
    5. Angle at which both rhoptries are docked to the AV, compared to the long axis of the AV
    
**This script differs from AV_measurements.py in that it utilizes pandas dataframes and only requires input of model points from model2point -object 
    
USAGE:
    Run the script in the directory where all your *AV.mod files are located
"""

import os
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.stats
import re
import pandas as pd
import seaborn as sns
sns.set()

#Run a shell command to convert .mod files to .txt files of the same name 
os.system('ls *AV.mod | awk -F "." \'{print "model2point -object "$1".mod "$1"_ObjPts.txt"}\' | sh')

#Read all .txt files for models into list pointFiles
stream = os.popen('ls *ObjPts.txt').read()
pointFiles = stream.split("\n"); pointFiles.pop()

allModels = {}

#Read model points into a pandas dataframe, add all these to a dict (allModels) with tomo ID as key
for x in pointFiles:
    with open(x,'r') as f:
        points = f.readlines()
    points = [re.split(r'\s+',line.strip()) for line in points] #Seems like re.split could be substitutes with just line.strip().split()
    
    #Convert points into a Pandas dataframe
    df = pd.DataFrame({'Ob':[int(points[0][0])], 'Ct':[int(points[0][1])], 'X':[float(points[0][2])], 'Y':[float(points[0][3])], 'Z':[float(points[0][4])]})
    count = 1 
    for line in points[1:]:
        df = df.append(pd.DataFrame({'Ob':[int(points[count][0])], 'Ct':[int(points[count][1])], 'X':[float(points[count][2])], 'Y':[float(points[count][3])], 'Z':[float(points[count][4])]}), ignore_index=True)
        count += 1
        
    nameTomo = os.path.basename(x).replace('_AV_ObjPts.txt','') #replaces end of model file name with nothing, so nameTomo is the tomo ID
    allModels[nameTomo] = df

pixSize = 1.0616 #Pixel size in nanometers

#Recompile AV properties as follows:
AvProps = None
AvProps = pd.DataFrame(columns = ['Name', 'Dist', 'MajAxis', 'MinAxis', 'Ecc', 'AngleAV', 'AngleAV2', 'AngleAV3', 'AngleRh1', 'AngleRh2', 'NumRhop', 'RhopDist1', 'RhopDist2', 'RhopGuided'])

#Two functions to calculate distances between points
def calcDist(a,b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def calcTotalDist(a):
    cumDist = 0
    for i in range(len(a)-1):
        dist = calcDist(a[i],a[i+1])
        cumDist += dist
    return cumDist

#More info on the next two functions can be found in python file from Shrawan called ShortestDistance.py
def t(p,q,r):
    x = p-q
    return np.dot(r-q,x)/np.dot(x,x)

def calcShortest(p,q,r): #Points p and q define the line. Calculate the shortest distance from r to line p-q
    return np.linalg.norm(t(p,q,r)*(p-q)+q-r) #linalg.norm is the same as calcDist()

def calcShortestIntersection(p,q,r):
    return t(p,q,r)*(p-q)


def line_intersection(line1, line2): #Find intersection of 2 lines. Sidenote: I literally do not know the math that is happening here. Copied from Shrawan's python script Av2_toxo.py
    OB = line1[1]
    OA = line1[0]
    AB = (OB[0]-OA[0],OB[1]-OA[1],OB[2]-OA[2])
    # parametric equation of AB r = OA + t(AB)
    r1 = ((OA[0],AB[0]),(OA[1],AB[1]),(OA[2],AB[2]))

    OD = line2[1]
    OC = line2[0]
    CD = (OD[0]-OC[0],OD[1]-OC[1],OD[2]-OC[2])
    # parametric equation of CD r = OC + s(CD)
    r2 = ((OC[0],CD[0]),(OC[1],CD[1]),(OC[2],CD[2]))

    # r1=r2, bring the 3 equations to the form a*t + b*s = c and convert to matrix
    # form of Az = C 

    # 3 equations
    #r1[0][1] - r2[0][1] = r2[0][0]-r1[0][0]
    #r1[1][1] - r2[1][1] = r2[1][0]-r1[1][0]
    #r1[2][1] - r2[2][1] = r2[2][0]-r1[2][0]
    A = np.array([(r1[0][1],-1*r2[0][1]),(r1[1][1],-1*r2[1][1]),
                   (r1[2][1],-1*r2[2][1])])
    C = np.array([r2[0][0]-r1[0][0],r2[1][0]-r1[1][0],r2[2][0]-r1[2][0]])

    ts = np.linalg.lstsq(A,C) # returns values for (t,s) and sum of residuals etc
    #    print(ts)
    x_int = r1[0][0] + (r1[0][1]*ts[0][0])
    y_int = r1[1][0] + (r1[1][1]*ts[0][0])
    z_int = r1[2][0] + (r1[2][1]*ts[0][0])

    return np.array([x_int,y_int,z_int])

def calcAngle(a,b,c):
    ba = a-b
    bc = c-b
    cosAngle = np.dot(ba, bc)/(np.linalg.norm(ba)*np.linalg.norm(bc))
    angle = np.arccos(cosAngle)*180/np.pi
    return angle

#Allocate all points to different features: Apex, AV axes, Rhoptries
for key in allModels:
    modelPoints = allModels[key]
    
    #Apex = central density of RSA (basically the closts point between PPM and AV)
    apex = modelPoints.loc[modelPoints['Ob'] == 1] #Grabs rows where the object number is 1

    if len(apex) == 2:
        AvFront = apex[['X','Y','Z']].values[1] #1D array
        apex = apex[['X','Y','Z']].values[0] #1D array
    else:
        print('Apex of parasite %s is absent or improperly marked' %key)
        AvFront = np.array(None)
        apex = np.array(None)
        
    #Create a unique AV name for each key
    AvName = key + '_AV'
    
    #Major Axis of AV
    majorAxis = modelPoints.loc[(modelPoints['Ob'] == 2) & (modelPoints['Ct'] == 1)]
    if len(majorAxis) == 2:
        majorAxis = majorAxis[['X','Y','Z']].values #2D array
    else:
        print('Major axis of AV %s is absent or has improper contour' %AvName)
        majorAxis = np.array(None)
        
    #Minor Axis of AV
    minorAxis = modelPoints.loc[(modelPoints['Ob'] == 2) & (modelPoints['Ct'] == 2)]
    if len(minorAxis) == 2:
        minorAxis = minorAxis[['X','Y','Z']].values #2D array
    else:
        print('Minor axis of AV %s is absent or has improper contour' %AvName)
        minorAxis = np.array(None)
        
    #Point where rhoptry most parallel to AV is closest to the AV
    RhopCon1 = modelPoints.loc[(modelPoints['Ob'] == 3) & (modelPoints['Ct'] == 1)]
    if len(RhopCon1) == 1:
        apexRhop1 = RhopCon1[['X','Y','Z']].values[0] #1D array
    else:
        print('Rhoptry 1 of AV %s is absent or improperly marked' %AvName)
        apexRhop1 = np.array(None)
        
    #Point where second rhoptry docks to the AV
    RhopCon2 = modelPoints.loc[(modelPoints['Ob'] == 3) & (modelPoints['Ct'] == 2)]
    if len(RhopCon2) == 2:
        apexRhop2 = RhopCon2[['X','Y','Z']].values[0] #1D array
        AvRhop2 = RhopCon2[['X','Y','Z']].values[1] #1D array
    else:
        print('Rhoptry 2 of AV %s is absent of improperly marked' %AvName)
        apexRhop2 = np.array(None)
        AvRhop2 = np.array(None)
        
    #All calculations
    #Major Axis in nm
    if majorAxis.all() != None:
        majAxisLen = calcDist(majorAxis[0],majorAxis[1])*pixSize
    else:
        majAxisLen = np.nan
    
    #Minor axis length in nm
    if minorAxis.all() != None:
        minAxisLen = calcDist(minorAxis[0],minorAxis[1])*pixSize
    else:
        minAxisLen = np.nan
    
    #Swap coordinates for minor and major axis if the minor axis is longer than the major axis
    if majAxisLen < minAxisLen:
        minorAxis, majorAxis = majorAxis, minorAxis
        
        if majorAxis.all() != None:
            majAxisLen = calcDist(majorAxis[0],majorAxis[1])*pixSize
        else:
            majAxisLen = np.nan
        
        #Minor axis length in nm
        if minorAxis.all() != None:
            minAxisLen = calcDist(minorAxis[0],minorAxis[1])*pixSize
        else:
            minAxisLen = np.nan
    
    #Eccentricity of AV
    if (majAxisLen != np.nan) and (minAxisLen != np.nan):
        b = max(majAxisLen,minAxisLen)
        a = min(majAxisLen,minAxisLen)
        ecc = math.sqrt((b/2)**2-(a/2)**2)/(b/2)
        
    #Distance from apex (central density of RSA) to AV
    if (AvFront.all() != None) and (apex.all() != None):
        AvDist = calcDist(apex,AvFront)*pixSize
    else:
        AvDist = np.nan
        
    #Calculate shortest distance from apex to major axis line (extended)
    #Calculate angle at which AV major axis is oriented with respect to the apex, using the AV center point
    if (majorAxis.all() != None) and (apex.all() != None):
        shortestDist1 = calcShortest(majorAxis[1], majorAxis[0], apex)*pixSize
        AvCenter = line_intersection(majorAxis, minorAxis)
        hypot = calcDist(AvCenter, apex)*pixSize
        angleAV1 = np.arcsin(shortestDist1/hypot)*180/np.pi #Multipling by 180/pi converts to degrees
    else:
        angleAV1 = np.nan
    
    #Calculate shortest distance frmo apex to major axis line
    #Calculate angle at which the AV major axis is oriented with respect to the apex, using the AV-apex line
    if (majorAxis.all() != None) and (apex.all() != None):
        intersect = calcShortestIntersection(majorAxis[1], majorAxis[0], apex)
        angle = calcAngle(AvFront, apex, intersect)
        angleAV2 = 180-90-angle
    else:
        angleAV2 = np.nan
        
    #Calculate intersection between Apex-AvFront line and major axis line, then calculate the angle between the two
    if (majorAxis.all() != None) and (apex.all() != None):
        apexLine = modelPoints.loc[modelPoints['Ob'] == 1]
        apexLine = apexLine[['X','Y','Z']].values
        intersect2 = line_intersection(majorAxis, apexLine)
        hypot2 = calcDist(intersect2, apex)
        angleAV3 = np.arcsin(shortestDist1/hypot2)*180/np.pi #Calculate in degrees
    else:
        angleAV3 = np.nan
    
    #Calculate angle and distance of rhoptry 1 docking with respect to the AV major axis
    if (apexRhop1.all() != None) and (majorAxis.all() != None):
        angleRhop1 = calcAngle(apexRhop1, majorAxis[0], majorAxis[1])
        distRhop1 = calcDist(apexRhop1, majorAxis[1])*pixSize
    else:
        angleRhop1 = np.nan
        distRhop1 = np.nan
        
    #Calculate angle of rhoptry 2 docking with respect to the AV major axis
    if (apexRhop2.all() != None) and (majorAxis.all() != None):
        angleRhop2 = calcAngle(apexRhop2, majorAxis[0], majorAxis[1])
    else:
        angleRhop2 = np.nan
    
    #Calculate distance between rhoptry 2 and AV
    if (apexRhop2.all() != None) and (AvRhop2.all() != None):
        distRhop2 = calcDist(apexRhop2, AvRhop2)*pixSize
    else:
        distRhop2 = np.nan
    
    #Quantify the number of rhoptries docked at the AV
    RhopCount = 0
    """
    if distRhop1 <= 50:
        RhopCount += 1
    if distRhop2 <= 50:
        RhopCount += 1
    """
      
    RhopCount = 0
    if math.isnan(distRhop1) != True:
        RhopCount += 1
    if math.isnan(distRhop2) != True:
        RhopCount += 1
    
    
    #Quantify rhoptries guided to the apex, but not docked, within a certain distance range
    RhopGuid = 0
    if 50 < distRhop1 <= 100:
        RhopGuid += 1
    if 50 < distRhop2 <= 100:
        RhopGuid += 1
        
    if (RhopCount == 2) or (RhopCount == 1 and RhopGuid == 0):
        RhopGuid = float("nan")
    
    
    #Reclassify rhoptry 1 vs rhoptry 2 more accureately based on angle 
    for i in range(len(AvProps['AngleRh1'].values)):
        if math.isnan(AvProps['AngleRh1'].values[i]) == True or math.isnan(AvProps['AngleRh2'].values[i]) == True:
            continue
        elif AvProps['AngleRh2'].values[i] < AvProps['AngleRh1'].values[i]:
            AvProps['AngleRh1'].values[i], AvProps['AngleRh2'].values[i] = AvProps['AngleRh2'].values[i], AvProps['AngleRh1'].values[i] #Swaps the 2 values between lists if angle of rhoptry 2 is less than angle of rhoptry 1
            AvProps['RhopDist1'].values[i], AvProps['RhopDist2'].values[i] = AvProps['RhopDist2'].values[i], AvProps['RhopDist1'].values[i] #Also swaps the rhoptry distances if it swaps the angles
        else:
            continue
    '''
    
    #Reclassify rhoptry 2 vs rhoptry 2 based on distance from AV. This should not be used at the same time
    #as the above chunk of code, which classifies rhoptries 1 vs 2 based on angle
    for i in range(len(AvProps['RhopDist1'].values)):
        if math.isnan(AvProps['RhopDist1'].values[i]) == True or math.isnan(AvProps['RhopDist2'].values[i]) == True:
            continue
        elif AvProps['RhopDist2'].values[i] < AvProps['RhopDist1'].values[i]:
            AvProps['RhopDist1'].values[i], AvProps['RhopDist2'].values[i] = AvProps['RhopDist2'].values[i], AvProps['RhopDist1'].values[i] #Swaps the 2 values between lists if distance of rhoptry 2 is less than distance of rhoptry 1
            AvProps['AngleRh1'].values[i], AvProps['AngleRh2'].values[i] = AvProps['AngleRh2'].values[i], AvProps['AngleRh1'].values[i] #Also swaps rhoptry angles if it swaps the distances
        else:
            continue
    '''
        
    #Create a new dataframe comprising all the properties
    AvProps = AvProps.append(pd.DataFrame({'Name':[AvName], 'Dist':[AvDist], 'MajAxis':[majAxisLen], 'MinAxis':[minAxisLen], 'Ecc':[ecc], 'AngleAV':[angleAV1], 'AngleAV2':[angleAV2], 'AngleAV3':[angleAV3], 'AngleRh1':[angleRhop1], 'AngleRh2':[angleRhop2], 'NumRhop':[RhopCount], 'RhopDist1':[distRhop1], 'RhopDist2':[distRhop2], 'RhopGuided':[RhopGuid]}), ignore_index = True)
    
print ('Mean eccentricity of Av is', np.mean(AvProps['Ecc']), 'and a std of', np.std(AvProps['Ecc']))

#PLOT ALL THE DATA

Figure1 = plt.figure(1, figsize=(4,4))
ax1 = plt.axes([0.3,0.05,0.65,0.9])
ax1.set_ylabel('Av dist (nm)',fontsize=30, labelpad=8)
#ax1.set_axis_bgcolor(plt.cm.Blues(0.05))
ax1.grid(False)
ax1.tick_params('x', top=False, bottom=False, labelsize=0)
ax1.tick_params('y', left=True, labelsize=16)
ax1.set_ylim(-1,34)#max(AvProps['Dist'])+5)
ax1.set_xlim(-0.3,0.3)
meanAvDist = np.mean(AvProps['Dist'])
stdAvDist = np.std(AvProps['Dist'])
ax1.plot([-0.53,0.53],[meanAvDist,meanAvDist],color='k',ls='-',lw=5,alpha=0.7)
ax1 = sns.swarmplot(data=AvProps['Dist'],s=12,color='b',alpha=0.9)
markers, caps, bars = ax1.errorbar(x=0,y=meanAvDist,yerr=stdAvDist,elinewidth=5,
                                                                ecolor='r',
                                                                capsize=30,capthick=5)  

[bar.set_alpha(0.99) for bar in bars]
[cap.set_alpha(0.99) for cap in caps]   
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/AvDist.png')   

Figure2 = plt.figure(2, figsize=(4,4))
ax2 = plt.axes([0.25,0.25,0.70,0.7])
ax2.plot(AvProps['MinAxis'], AvProps['MajAxis'], marker='o', mfc='b', mec='k',
         markersize=12, lw=0, alpha=0.4) # swap x and y axes for toxo data
markers, caps, bars = ax2.errorbar(x=np.mean(AvProps['MinAxis']),y=np.mean(AvProps['MajAxis']),
                       xerr=np.std(AvProps['MinAxis']),yerr=np.std(AvProps['MajAxis']),
                                   ecolor='r',elinewidth=4,capsize=5,capthick=4)
[bar.set_alpha(0.99) for bar in bars]
[cap.set_alpha(0.99) for cap in caps]
#markers.set_alpha(0.99) 

ax2.set_xlim(-2,max(max(AvProps['MinAxis']), max(AvProps['MajAxis']))+5)
ax2.set_ylim(-2,max(max(AvProps['MinAxis']), max(AvProps['MajAxis']))+5)
ax2.tick_params(labelsize=16)
ax2.set_xlabel('Av min (nm)',fontsize=30, labelpad=8)
ax2.set_ylabel('Av maj (nm)',fontsize=30, labelpad=6)
#Plot line with slope of 1
x_vals = np.array(ax2.get_xlim())
y_vals = 1*x_vals+0
ax2.plot(x_vals,y_vals,'--')
#ax2.set_axis_bgcolor(plt.cm.Blues(0.05))
ax2.grid(False)
ax2.tick_params('x', bottom=True, labelsize=16)
ax2.tick_params('y', left=True)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/AvDimensions.png')


Figure3 = plt.figure(3, figsize=(4,4))
ax3 = plt.axes([0.3,0.05,0.65,0.9])
ax3.set_ylabel(r'$\psi$',fontsize=30, labelpad=8)
#ax3.set_axis_bgcolor(plt.cm.Blues(0.05))
ax3.grid(False)
ax3.tick_params('x', top=False, bottom=False, labelsize=0)
ax3.tick_params('y', right=False, labelsize=16)
ax3.set_ylim(-1,max(AvProps['AngleAV'])+5)
ax3.set_xlim(-0.3,0.3)

#Remove NaN values
noNaNAv = AvProps.loc[np.logical_not(np.isnan(AvProps['AngleAV']))]['AngleAV']
medAngleAv = np.median(noNaNAv)
#ax3.plot([-0.53,0.53],[medAngleAv,medAngleAv],color='k',ls='-',lw=5,alpha=0.7)

ax3 = sns.swarmplot(data=AvProps['AngleAV'],s=12,color='b',alpha=0.9)
sns.boxplot(data=AvProps['AngleAV'], color='w', ax=ax3)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/AvAngle.png')

'''
Figure4 = plt.figure(4, figsize=(4,4))
ax4 = plt.axes([0.3,0.05,0.65,0.9])
ax4.set_ylabel(r'$\psi$',fontsize=30, labelpad=8)
#ax3.set_axis_bgcolor(plt.cm.Blues(0.05))
ax4.grid(False)
ax4.tick_params('x', top=False, bottom=False, labelsize=0)
ax4.tick_params('y', right=False, labelsize=16)
ax4.set_ylim(-1,max(AvProps['AngleAV2'])+5)
ax4.set_xlim(-0.3,0.3)

#Remove NaN values
noNaNAv2 = AvProps.loc[np.logical_not(np.isnan(AvProps['AngleAV2']))]['AngleAV2']
medAngleAv2 = np.median(noNaNAv2)
ax4.plot([-0.53,0.53],[medAngleAv2,medAngleAv2],color='k',ls='-',lw=5,alpha=0.7)
ax4 = sns.swarmplot(data=AvProps['AngleAV2'],s=12,color='b',alpha=0.9)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/AvAngle2.png')


Figure5 = plt.figure(5, figsize=(4,4))
ax5 = plt.axes([0.3,0.05,0.65,0.9])
ax5.set_ylabel(r'$\psi$',fontsize=30, labelpad=8)
#ax3.set_axis_bgcolor(plt.cm.Blues(0.05))
ax5.grid(False)
ax5.tick_params('x', top=False, bottom=False, labelsize=0)
ax5.tick_params('y', right=False, labelsize=16)
ax5.set_ylim(-1,max(AvProps['AngleAV3'])+5)
ax5.set_xlim(-0.3,0.3)

#Remove NaN values
noNaNAv3 = AvProps.loc[np.logical_not(np.isnan(AvProps['AngleAV3']))]['AngleAV3']
medAngleAv3 = np.median(noNaNAv3)
ax5.plot([-0.53,0.53],[medAngleAv3,medAngleAv3],color='k',ls='-',lw=5,alpha=0.7)
ax5 = sns.swarmplot(data=AvProps['AngleAV3'],s=12,color='b',alpha=0.9)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/AvAngle3.png')
'''

Figure6 = plt.figure(6, figsize=(4,4))
ax6 = plt.axes([0.3,0.05,0.65,0.9])
ax6.set_ylabel(r'$\theta_{Rh1}$',fontsize=30, labelpad=8)
#ax6.set_axis_bgcolor(plt.cm.Blues(0.05))
ax6.grid(False)
ax6.tick_params('x', top=False, bottom=False, labelsize=0)
ax6.tick_params('y', right=False, labelsize=16)
ax6 = sns.swarmplot(data=AvProps['AngleRh1'],s=12,color='b',alpha=0.9)

#NaN values removed
noNaNRh1 = AvProps.loc[np.logical_not(np.isnan(AvProps['AngleRh1']))]['AngleRh1']
medAngleRh1 = np.median(noNaNRh1)
#ax6.plot([-0.53,0.53],[medAngleRh1,medAngleRh1],color='k',ls='-',lw=5,alpha=0.7)
ax6.set_ylim(-1,max(AvProps['AngleRh1'].values)+5)
ax6.set_xlim(-0.3,0.3)
sns.boxplot(data=AvProps['AngleRh1'], color='w',ax=ax6)

#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/RhAngle1.png')


Figure7 = plt.figure(7, figsize=(4,4))
ax7 = plt.axes([0.3,0.05,0.65,0.9])
ax7.set_ylabel(r'$\theta_{Rh2}$',fontsize=30, labelpad=8)
#ax7.set_axis_bgcolor(plt.cm.Blues(0.05))
ax7.grid(False)
ax7.tick_params('x', top=False, bottom=False, labelsize=0)
ax7.tick_params('y', right=False, labelsize=16)
ax7 = sns.swarmplot(data=AvProps['AngleRh2'],s=12,color='b',alpha=0.9)

#NaN values removed
noNaNRh2 = AvProps.loc[np.logical_not(np.isnan(AvProps['AngleRh2']))]['AngleRh2']
medAngleRh2 = np.median(noNaNRh2)
#ax7.plot([-0.53,0.53],[medAngleRh2,medAngleRh2],color='k',ls='-',lw=5,alpha=0.7)
ax7.set_ylim(-1,max(AvProps['AngleRh2'].values)+5)
ax7.set_xlim(-0.3,0.3)
sns.boxplot(data=AvProps['AngleRh2'], color='w',ax=ax7)

#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/RhAngle2.png')



Figure8 = plt.figure(8, figsize=(4,4))
ax8 = plt.axes([0.23,0.35,0.72,0.60])
rhop0 = len(AvProps.loc[AvProps['NumRhop']==0])
rhop1 = len(AvProps.loc[AvProps['NumRhop']==1])
rhop2 = len(AvProps.loc[AvProps['NumRhop']==2])
ax8.bar([0,1,2],[rhop0,rhop1,rhop2],align='center',tick_label=['0','1','2'],
        color=('b'))
ax8.set_xlim(-0.7,2.7)
ax8.set_ylim(0,max(rhop0,rhop1,rhop2)+2)
ax8.tick_params(labelsize=16)
ax8.set_xlabel('Rhoptries docked \non Av',fontsize=30, labelpad=8)
ax8.set_ylabel('Frequency (#)',fontsize=30, labelpad=6)
#ax8.set_axis_bgcolor(plt.cm.Blues(0.05))
ax8.grid(False)
ax8.tick_params('x', top=False, labelsize=16)
ax8.tick_params('y', right=False)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/RhDocked.png')


Figure9 = plt.figure(9, figsize=(4,4))
ax9 = plt.axes([0.25,0.25,0.70,0.7])
ax9.plot(AvProps['AngleRh1'], AvProps['AngleRh2'], marker='o', mfc='b', mec='k',
         markersize=12, lw=0, alpha=0.4) 
markers, caps, bars = ax9.errorbar(x=np.mean(noNaNRh1),y=np.mean(noNaNRh2),
                       xerr=np.std(noNaNRh1),yerr=np.std(noNaNRh2),
                                   ecolor='r',elinewidth=4,capsize=5,capthick=4)
[bar.set_alpha(0.99) for bar in bars]
[cap.set_alpha(0.99) for cap in caps]
#markers.set_alpha(0.99) 


ax9.set_xlim(-2,max(max(AvProps['AngleRh2'].values), max(AvProps['AngleRh1'].values))+5)
ax9.set_ylim(-2,max(max(AvProps['AngleRh2'].values), max(AvProps['AngleRh1']))+5)
ax9.tick_params(labelsize=16)
ax9.set_xlabel(r'$\theta_{Rh1}$',fontsize=30, labelpad=8)
ax9.set_ylabel(r'$\theta_{Rh2}$',fontsize=30, labelpad=6)

#Plot line with slope of 1
x_vals = np.array(ax9.get_xlim())
y_vals = 1*x_vals+0
ax9.plot(x_vals,y_vals,'--')
#ax9.set_axis_bgcolor(plt.cm.Blues(0.05))
ax9.grid(False)
ax9.tick_params('x', bottom=True, labelsize=16)
ax9.tick_params('y', left=True)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/RhAngleVsAngle.png')


Figure10 = plt.figure(10, figsize=(4,4))
ax10 = plt.axes([0.25,0.25,0.70,0.7])
ax10.plot(AvProps['AngleAV'], AvProps['Dist'], 'o', mfc='b', mec='k', markersize=12,
         alpha=0.4)

#NaN values removed
noNaNDist1 = AvProps.loc[np.logical_not(np.isnan(AvProps['Dist']))]['Dist']
#Plot best fit line
m, b = np.polyfit(noNaNAv, noNaNDist1, 1)
ax10.plot(noNaNAv, m*(noNaNAv)+b)
slope1, intercept1, r1, p1, se1 = scipy.stats.linregress(noNaNAv, noNaNDist1)
ax10.set_xlim(-2,max(AvProps['AngleAV'])+5)
ax10.set_ylim(-2,max(AvProps['Dist'])+5)
ax10.tick_params(labelsize=16)
ax10.set_xlabel(r'$\psi$',fontsize=30, labelpad=8)
ax10.set_ylabel('Av dist (nm)',fontsize=30, labelpad=6)
#ax10.set_axis_bgcolor(plt.cm.Blues(0.05))
ax10.grid(False)
ax10.tick_params('x', top=False, labelsize=16)
ax10.tick_params('y', right=False)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/AvAngle1VsDist.png')


Figure11 = plt.figure(11, figsize=(4,4))
ax11 = plt.axes([0.3,0.05,0.65,0.9])
ax11.set_ylabel('Rhoptry1 dist (nm)',fontsize=30, labelpad=8)
#ax1.set_axis_bgcolor(plt.cm.Blues(0.05))
ax11.grid(False)
ax11.tick_params('x', top=False, bottom=False, labelsize=0)
ax11.tick_params('y', left=True, labelsize=16)
ax11.set_ylim(-1,max(AvProps['RhopDist1'])+5)
ax11.set_xlim(-0.3,0.3)
meanRhopDist1 = np.mean(AvProps['RhopDist1'])
stdRhopDist1 = np.std(AvProps['RhopDist1'])
ax11.plot([-0.53,0.53],[meanRhopDist1,meanRhopDist1],color='k',ls='-',lw=5,alpha=0.7)
ax11 = sns.swarmplot(data=AvProps['RhopDist1'],s=12,color='b',alpha=0.9)
markers, caps, bars = ax11.errorbar(x=0,y=meanRhopDist1,yerr=stdRhopDist1,elinewidth=5,
                                                                ecolor='r',
                                                                capsize=30,capthick=5)
[bar.set_alpha(0.99) for bar in bars]
[cap.set_alpha(0.99) for cap in caps]  
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/RhDist1.png')


Figure12 = plt.figure(12, figsize=(4,4))
ax12 = plt.axes([0.25,0.25,0.70,0.7])
ax12.plot(AvProps['AngleRh1'], AvProps['RhopDist1'], 'o', mfc='b', mec='k', markersize=12,
         alpha=0.4)

#NaN values removed
noNaNRhDist1 = AvProps.loc[np.logical_not(np.isnan(AvProps['RhopDist1']))]['RhopDist1']
#Plot best fit line
m, b = np.polyfit(noNaNRh1, noNaNRhDist1, 1)
ax12.plot(noNaNRh1, m*(noNaNRh1)+b)
slope2, intercept2, r2, p2, se2 = scipy.stats.linregress(noNaNRh1, noNaNRhDist1)
ax12.set_xlim(-2,max(AvProps['AngleRh1'])+5)
ax12.set_ylim(-2,max(AvProps['RhopDist1'])+5)
ax12.tick_params(labelsize=16)
ax12.set_xlabel(r'$\theta_{Rh1}$',fontsize=30, labelpad=8)
ax12.set_ylabel('Rhoptry1 dist (nm)',fontsize=30, labelpad=6)
#ax10.set_axis_bgcolor(plt.cm.Blues(0.05))
ax12.grid(False)
ax12.tick_params('x', top=False, labelsize=16)
ax12.tick_params('y', right=False)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/RhAngle1vsRhDist1.png')


Figure13 = plt.figure(13, figsize=(4,4))
ax13 = plt.axes([0.3,0.05,0.65,0.9])
ax13.set_ylabel('Rhoptry2 dist (nm)',fontsize=30, labelpad=8)
#ax1.set_axis_bgcolor(plt.cm.Blues(0.05))
ax13.grid(False)
ax13.tick_params('x', top=False, bottom=False, labelsize=0)
ax13.tick_params('y', left=True, labelsize=16)
ax13.set_ylim(-1,max(AvProps['RhopDist2'])+5)
ax13.set_xlim(-0.3,0.3)
meanRhopDist2 = np.mean(AvProps['RhopDist2'])
stdRhopDist2 = np.std(AvProps['RhopDist2'])
ax13.plot([-0.53,0.53],[meanRhopDist2,meanRhopDist2],color='k',ls='-',lw=5,alpha=0.7)
ax13 = sns.swarmplot(data=AvProps['RhopDist2'],s=12,color='b',alpha=0.9)
markers, caps, bars = ax13.errorbar(x=0,y=meanRhopDist2,yerr=stdRhopDist2,elinewidth=5,
                                                                ecolor='r',
                                                                capsize=30,capthick=5)
[bar.set_alpha(0.99) for bar in bars]
[cap.set_alpha(0.99) for cap in caps]  
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/RhDist2.png')


Figure14 = plt.figure(14, figsize=(4,4))
ax14 = plt.axes([0.25,0.25,0.70,0.7])
ax14.plot(AvProps['AngleRh2'], AvProps['RhopDist2'], 'o', mfc='b', mec='k', markersize=12,
         alpha=0.4)

#NaN values removed
noNaNRhDist2 = AvProps.loc[np.logical_not(np.isnan(AvProps['RhopDist2']))]['RhopDist2']
#Plot best fit line
m, b = np.polyfit(noNaNRh2, noNaNRhDist2, 1)
ax14.plot(noNaNRh2, m*(noNaNRh2)+b)
slope3, intercept3, r3, p3, se3 = scipy.stats.linregress(noNaNRh2, noNaNRhDist2)
ax14.set_xlim(-2,max(AvProps['AngleRh2'])+5)
ax14.set_ylim(-2,max(AvProps['RhopDist2'])+5)
ax14.tick_params(labelsize=16)
ax14.set_xlabel(r'$\theta_{Rh2}$',fontsize=30, labelpad=8)
ax14.set_ylabel('Rhoptry2 dist (nm)',fontsize=30, labelpad=6)
#ax10.set_axis_bgcolor(plt.cm.Blues(0.05))
ax14.grid(False)
ax14.tick_params('x', top=False, labelsize=16)
ax14.tick_params('y', right=False)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/RhAngle2vsRhDist2.png')


Figure15 = plt.figure(15, figsize=(4,4))
ax15 = plt.axes([0.25,0.25,0.70,0.7])
ax15.plot(AvProps['RhopDist1'], AvProps['RhopDist2'], marker='o', mfc='b', mec='k',
         markersize=5, lw=0)#, alpha=0.4) 
#markers, caps, bars = ax15.errorbar(x=np.mean(noNaNRhDist1),y=np.mean(noNaNRhDist2),
                       #xerr=np.std(noNaNRhDist1),yerr=np.std(noNaNRhDist2),
                                   #ecolor='r',elinewidth=2,capsize=5,capthick=2)
#[bar.set_alpha(0.99) for bar in bars]
#[cap.set_alpha(0.99) for cap in caps]
#markers.set_alpha(0.99)


ax15.set_xlim(-5,120)#(-2,max(max(AvProps['RhopDist1'].values), max(AvProps['RhopDist2'].values))+5)
ax15.set_ylim(-5,120)#(-2,max(max(AvProps['RhopDist1'].values), max(AvProps['RhopDist2']))+5)
ax15.tick_params(labelsize=16)
ax15.set_xlabel('Rhoptry1 dist (nm)',fontsize=30, labelpad=8)
ax15.set_ylabel('Rhoptry2 dist (nm)',fontsize=30, labelpad=6)

x_vals = np.array(ax15.get_xlim())
y_vals = 1*x_vals+0
#ax15.plot(x_vals,y_vals,'--')
#ax15.set_axis_bgcolor(plt.cm.Blues(0.05))
ax15.grid(False)
ax15.tick_params('x', bottom=True,)
ax15.tick_params('y', left=True)
plt.axvline(50,ymax = 1, color = 'r', linestyle = '--')
plt.axvline(100,ymax = 1, color = 'r', linestyle = '--')
plt.axhline(50,xmax = 1, color = 'r', linestyle = '--')
plt.axhline(100,xmax = 1, color = 'r', linestyle = '--')

ax15.fill_between([-5,50],50,-5,color='C6',alpha=0.3)
ax15.fill_between([50,100],50,-5,color='C1',alpha=0.3)
ax15.fill_between([-5,50],100,50,color='C1',alpha=0.3)
ax15.fill_between([50,100],100,50,color='C2',alpha=0.3)
ax15.fill_between([-5,50],120,100,color='C7',alpha=0.3)
ax15.fill_between([100,120],50,-5,color='C7',alpha=0.3)
ax15.fill_between([50,100],120,100,color='C8',alpha=0.3)
ax15.fill_between([100,120],100,50,color='C8',alpha=0.3)

#ax15.fill_between([-5,120],120,-5,color='C3',alpha=0.3,zorder=0)
#ax15.fill_between([100,120],120,-5,color='C3',alpha=0.3)

#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Plasmodium/ApicalVesicle/Analysis/Figures/RhDist1vsRhDist2.png')


Figure16 = plt.figure(16, figsize=(4,4))
ax16 = plt.axes([0.23,0.35,0.72,0.60])
rhopGuide0 = len(AvProps.loc[AvProps['RhopGuided']==0])
rhopGuide1 = len(AvProps.loc[AvProps['RhopGuided']==1])
rhopGuide2 = len(AvProps.loc[AvProps['RhopGuided']==2])

ax16.bar([0,1,2],[rhopGuide0,rhopGuide1,rhopGuide2],align='center',tick_label=['0','1','2'],
        color=('b'))
ax16.set_xlim(-0.7,2.7)
ax16.set_ylim(0,max(rhopGuide0,rhopGuide1,rhopGuide2)+2)
ax16.tick_params(labelsize=16)
ax16.set_xlabel('Rhoptries guided \nto apex',fontsize=30, labelpad=8)
ax16.set_ylabel('Frequency (#)',fontsize=30, labelpad=6)
#ax8.set_axis_bgcolor(plt.cm.Blues(0.05))
ax16.grid(False)
ax16.tick_params('x', top=False, labelsize=16)
ax16.tick_params('y', right=False)

Figure17 = plt.figure(17, figsize=(4,4))
ax17 = plt.axes([0.25,0.25,0.70,0.7])
ax17.plot(AvProps['RhopDist1'], AvProps['RhopDist2'], marker='o', mfc='b', mec='k', markersize=5, lw=0)
ax17.set_xlim(-5,120)#(-2,max(max(AvProps['RhopDist1'].values), max(AvProps['RhopDist2'].values))+5)
ax17.set_ylim(-5,120)#(-2,max(max(AvProps['RhopDist1'].values), max(AvProps['RhopDist2']))+5)
ax17.grid(False)
ax17.set_xlabel('Rhoptry1 dist (nm)',fontsize=30, labelpad=8)
ax17.set_ylabel('Rhoptry2 dist (nm)',fontsize=30, labelpad=6)
ax17.tick_params('x', bottom=True,)
ax17.tick_params('y', left=True)
plt.axvline(50,ymax = 1, color = 'r', linestyle = '--')
plt.axvline(100,ymax = 1, color = 'r', linestyle = '--')
plt.axhline(50,xmax = 1, color = 'r', linestyle = '--')
plt.axhline(100,xmax = 1, color = 'r', linestyle = '--')

ax17.fill_between([-5,50],50,-5,color='C6',alpha=0.3)
ax17.fill_between([50,100],50,-5,color='C1',alpha=0.3)
ax17.fill_between([-5,50],100,50,color='C1',alpha=0.3)
ax17.fill_between([50,100],100,50,color='C2',alpha=0.3)
ax17.fill_between([-5,50],120,100,color='C7',alpha=0.3)
ax17.fill_between([100,120],50,-5,color='C7',alpha=0.3)
ax17.fill_between([50,100],120,100,color='C8',alpha=0.3)
ax17.fill_between([100,120],100,50,color='C8',alpha=0.3)


#Definition for marginal axes
left, width = 0.25, 0.7
bottom, height = 0.25, 0.7
spacing = 0.05
rect_x = [left, bottom + height + spacing, width, 0.2]
rect_y = left + width + spacing, bottom, 0.2, height

ax17_x = Figure17.add_axes(rect_x, sharex=ax17)
ax17_x.grid(False)
ax17_x.axes.xaxis.set_visible(False)
ax17_x.axes.yaxis.set_visible(False)

ax17_y = Figure17.add_axes(rect_y, sharey=ax17)
ax17_y.axes.xaxis.set_visible(False)
ax17_y.axes.yaxis.set_visible(False)
ax17_y.grid(False)

sns.kdeplot(data=AvProps, y = 'RhopDist2', ax=ax17_y)
x_right = ax17_y.lines[-1].get_xdata()
y_right = ax17_y.lines[-1].get_ydata()
ax17_y.set_xlim(left=0,right=max(x_right)*1.2)
ax17_y.axhline(50,0,0.12, linestyle='-', color='r')
ax17_y.fill_betweenx(y_right,min(x_right),x_right,where=y_right<50,color='C6',alpha=0.3)
ax17_y.fill_betweenx(y_right,min(x_right),x_right,where=y_right>50,color='C1',alpha=0.3)

#ax17_y.fill_between([-5,50],  color='C6', alpha=0.3)
sns.kdeplot(data=AvProps, x = 'RhopDist1', ax=ax17_x)
x_top = ax17_x.lines[-1].get_xdata()
y_top = ax17_x.lines[-1].get_ydata()
ax17_x.set_ylim(bottom=0,top=max(y_top)*1.2)
ax17_x.axvline(50, 0, 0.04, linestyle='-', color='r')
ax17_x.fill_between(x_top,0,y_top,where=x_top < 50,color='C6',alpha=0.3)



