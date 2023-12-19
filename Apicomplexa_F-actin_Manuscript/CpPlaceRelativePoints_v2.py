#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 14:36:20 2022

@author: mattmartinez

Place points referring to a structure around another structure of interest.
This script will generate a ChimeraX command file. 

In this script for example, it will place points that refer to the tips of
actin filaments around a density of interest when opened in ChimeraX. Models 
are converted to .txt files before hand. This script acts on models with contours
that contain two (2) points.

It is important to note that the points were placed in tomograms, where the
structure of interest could be in a random rotation. I am plotting the data
on the subtomogram average, with a specific orientation. I will import alignment
data from the Dynamo refined table to rotate the points placed in the tomogram
"""

import numpy as np
from scipy.spatial.transform import Rotation as R
import scipy.stats as st
import glob
import re
import math
import subprocess
import csv
import mpl_toolkits.mplot3d as m3d
from matplotlib import cm
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

def calcDist(a,b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

Marker = [25,25,25] #Center point of the central PCR subunit
Marker2 = [float(29.6),float(14.9),float(14.6)] #Marker for the PD of the Crypto PCR. Point identified on Crypto PCR average in ChimeraX
Marker_TgPCR = [float(27.1),float(16.5),float(12.3)] #Marker for T. gondii knee density, from point on Toxo PCR average in ChimeraX
MarkerX = [float(44.5), float(12)] #2 X coordinates for correcting positions onto the central subunit. Used later in the script
PixSize = 10.6

allPoints = [] #Initialize list of points from the model files

model_folder            = "Models/"
file_list               = glob.glob("{}/*.txt".format(model_folder))
tbl                     = glob.glob("{}/*.tbl".format(model_folder))
tomoID_list = []
for item in file_list:
    item = item.split('_')
    tomoID_list.append(int(item[1]))

#allPoints is the points from the model files
for x in file_list:
    with open (x,'r') as f:
        points=f.readlines()
    points = [re.split(r'\s+', line.strip()) for line in points]
    allPoints.append(points)

#particles is the particles from the Dynamo refined table
with open (tbl[0],'r') as f:
    particles = f.readlines()
particles = [re.split(r'\s+', line.strip()) for line in particles]

particle_pairs = [] #First is particle from Dynamo, second is point from tomo
                    #List of 3: k (particle from Dynamo), i and counter (i + 
                    #counter is particle from model file)
PD_points = [] #List of points corresponding to the protruding density from
               #the tomogram models
actin_tomo_points = [] #List of points corresponding to the actin tip from the
                       #tomogram models

#Make lists of PD points, points for the actin tips
min_dist_list = []
for i in range(len(allPoints)):
    counter = 0
    while counter < len(allPoints[i]):
        min_dist = 100000
        point = [float(allPoints[i][counter][0]),float(allPoints[i][counter][1]),float(allPoints[i][counter][2])]
        PD_points.append(point)
        
        act_point = [float(allPoints[i][counter+1][0]),float(allPoints[i][counter+1][1]),float(allPoints[i][counter+1][2])]
        actin_tomo_points.append(act_point)
        
        tomo_particle_list = []
        for k in range(len(particles)):
            if tomoID_list[i] == int(particles[k][19]):
                tomo_particle_list.append(particles[k])
        
        for k in range(len(tomo_particle_list)):
            center = [float(tomo_particle_list[k][23])+float(tomo_particle_list[k][3]),float(tomo_particle_list[k][24])+float(tomo_particle_list[k][4]),float(tomo_particle_list[k][25])+float(tomo_particle_list[k][5])]
            dist = calcDist(center, act_point)
            if dist < min_dist:
                min_dist = dist
                pair = [int(tomo_particle_list[k][0]),i,counter]
            else:
                continue
        #print("particle {} with model {} particle {}".format(pair[0],pair[1],pair[2]))
        
        min_dist_list.append(min_dist)
        particle_pairs.append(pair)
        counter += 2
        
alignedPoints = []
centerPoints = []
EulerList = []
EulerFile = "{}/EulerAngles.csv".format(model_folder)
with open(EulerFile,'w', newline ='') as f:
    for i in range(len(particle_pairs)):
        center = [float(particles[particle_pairs[i][0]-1][23]),float(particles[particle_pairs[i][0]-1][24]),float(particles[particle_pairs[i][0]-1][25])]
        shift = [float(particles[particle_pairs[i][0]-1][3]),float(particles[particle_pairs[i][0]-1][4]),float(particles[particle_pairs[i][0]-1][5])]
        point = [(center[0]+shift[0]), (center[1]+shift[1]), (center[2]+shift[2])]
        alignedPoints.append(point)
        centerPoints.append(center)
        
        Eulers = [-float(particles[particle_pairs[i][0]-1][6]),-float(particles[particle_pairs[i][0]-1][7]),-float(particles[particle_pairs[i][0]-1][8])]
        EulerList.append(Eulers)
        write = csv.writer(f)
        write.writerow(Eulers)
        
    
#subprocess.run(["MOTL2Slicer", "{}/EulerAngles.csv".format(model_folder), "{}/SlicerAngles.csv".format(model_folder)])
SlicerFile = "{}/SlicerAngles.csv".format(model_folder)
with open(SlicerFile,'r') as f:
    angles = f.readlines()
angles = [re.split(r'\s+', line.strip()) for line in angles]

#Now, all refined center points for top PCR subunits and slicer angles are
#calculated and stores in lists. 
actin_vector_list = [] #A list for all vectors from aligned particle position 
                       #to the model point corresponding to the protruding density
                       #of the PCR. The aligned particle poisiton becomes the origin
                       #around which to rotate
Actin_points = []
Actin_points_onTgPCR = []
PD_positions = []
distances = []
actinKnee_distances = []
kneeCenter_distances = []

#Calculate vector and rotation to place actin tip in 3D
for i in range(len(alignedPoints)):
    angles[i] = [line.split(',') for line in angles[i]]
    angles[i] = [float(x) for x in angles[i][0]]
    
    shift = [float(particles[particle_pairs[i][0]-1][3]),float(particles[particle_pairs[i][0]-1][4]),float(particles[particle_pairs[i][0]-1][5])]
    vector = np.array([(actin_tomo_points[i][0]+shift[0])-alignedPoints[i][0],(actin_tomo_points[i][1]+shift[1])-alignedPoints[i][1],(actin_tomo_points[i][2]+shift[2])-alignedPoints[i][2]])
    
    #Rotate vector
    r = R.from_euler('ZXZ',[EulerList[i][0], EulerList[i][1], EulerList[i][2]], degrees=True)                       
    rotated_vector = r.apply(vector, inverse=True) 
    rotated_position = Marker + rotated_vector #Vector from the center point of the central PCR subunit
    
    dist = calcDist(rotated_position, Marker2)*1.06 #Distance from the rotated position to the PD (in nm)
    distances.append(dist)
    
    #Correct point if it is shifted onto neighboring subunit
    if abs(rotated_position[0]-MarkerX[0]) <= abs(rotated_position[0]-Marker2[0]):
        rotated_position = (rotated_position-[MarkerX[0], Marker2[1], Marker2[2]])+Marker2
    if abs(rotated_position[0]-MarkerX[1]) <= abs(rotated_position[0]-Marker2[0]):
        rotated_position = (rotated_position-[MarkerX[1], Marker2[1], Marker2[2]])+Marker2
        
    Actin_points.append(rotated_position)
    
    #To place Crypto points on the Toxo PCR subtomogram average
    vector_forTgPCR = (rotated_position - Marker2) + Marker_TgPCR
    Actin_points_onTgPCR.append(vector_forTgPCR)
    
    print("Vector: {}\nAngles: {}\nEuler Angles: {}\nRotated vector: {}\nActin point: {}\nRotated point: {}\n".format(vector, angles[i], EulerList[i], rotated_vector,act_point,rotated_position))

#From analagous script for Toxo, TgPlaceRelativePoints.py
TgDistances = [19.390044494547467, 39.09581951958778, 18.331401924390196, 23.812180559022032, 24.75514732302653, 38.94441269329018, 25.870537243534883, 40.499368064502434, 19.07162979131351, 45.93336794823061, 25.511421505807654, 41.40250653469919, 14.54015212543676, 28.58731488722518, 28.799743134535724, 44.000633384284804]


#Create command file that just places a marker at each point
with open('CpPlacePoints.cxc','w') as f:
    for i in range(len(Actin_points)):
        line = "marker #11 position {},{},{} color red radius 1.5\n".format(Actin_points[i][0],Actin_points[i][1],Actin_points[i][2])
        #line2 = "marker #11 position {},{},{} color blue radius 2\n".format(PD_positions[i][0],PD_positions[i][1],PD_positions[i][2])
        f.write(line)
        #f.write(line2)
        
with open('CpPlacePointsOnToxoPCR.cxc','w') as f:
    for i in range(len(Actin_points)):
        line = "marker #15 position {},{},{} color red radius 1.5\n".format(Actin_points_onTgPCR[i][0],Actin_points_onTgPCR[i][1],Actin_points_onTgPCR[i][2])
        f.write(line)

Cp_CI = st.t.interval(alpha=0.95, df=len(distances)-1, loc=np.mean(distances), scale=st.sem(distances))
Tg_CI = st.t.interval(alpha=0.95, df=len(TgDistances)-1, loc=np.mean(TgDistances), scale=st.sem(TgDistances))

Cp_Median = np.median(distances)
Tg_Median = np.median(TgDistances)

Cp_p75, Cp_p25 = np.percentile(distances, [75, 25])
Tg_p75, Tg_p25 = np.percentile(TgDistances, [75, 25])

distStatistic, distPvalue = st.mannwhitneyu(distances, TgDistances, alternative='two-sided')

#Plot of crypto distances
Figure1 = plt.figure(1, figsize=(4,4))
ax1 = plt.axes([0.3,0.05,0.5,0.9])
ax1.set_ylabel('Distance knee-actin tip (nm)', fontsize=20, labelpad=8)
ax1.set_facecolor(cm.Blues(0.05))
ax1.grid(False)
ax1.tick_params('x', top=False, bottom=False, labelsize=0)
ax1.tick_params('y', right=False, labelsize=16)
ax1.set_ylim(-5,max(distances)+10)

ax1_1 = sns.violinplot(data=distances, color=cm.Blues(0.2), width=0.5)
ax1_2 = sns.swarmplot(data=distances, s=8 ,color='b', alpha=0.9)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Toxoplasma/PCR/ActinKnee/Plots/CpKneeActinDist.png', dpi=600, bbox_inches='tight')

#Plot of crypto vs toxo distances
Figure2 = plt.figure(2, figsize=(4,4))
ax2 = plt.axes([0.3,0.05,0.5,0.9])
ax2.set_ylabel('Distance knee-actin tip (nm)', fontsize=20, labelpad=8)
ax2.set_facecolor(cm.Blues(0.05))
ax2.grid(False)
ax2.tick_params('x', top=False, bottom=False, labelsize=0)
ax2.tick_params('y', right=False, labelsize=16)
ax2.set_ylim(-5,max(distances)+10)

ax1_1 = sns.violinplot(data=[TgDistances,distances], color=cm.Blues(0.2), width=0.8)
ax1_2 = sns.swarmplot(data=[TgDistances,distances], s=6 ,color='b', alpha=0.9)
#plt.savefig('/Users/matthewmartinez/Desktop/Upenn/Yi-Wei_Chang_Lab/Toxoplasma/PCR/ActinKnee/Plots/CpVsTgKneeActinDist.png', dpi=600, bbox_inches='tight')
        
        
    
    


