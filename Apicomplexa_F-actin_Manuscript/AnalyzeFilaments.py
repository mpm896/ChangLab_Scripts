#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 13:29:16 2022

@author: mattmartinez
Plot interfilament angles and distances from the computational toolbax for filament analysis in MatLab
"""

import numpy as np
import scipy.stats as st
import random
import copy
import glob
import re
import math
import subprocess
import csv
from itertools import chain
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import seaborn as sns
sns.set()

my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap.colors[0])

#Import numbers from files
file_list = glob.glob('*.csv')

with open(file_list[0], newline = '') as f:
    reader = csv.reader(f)
    Jasp_Angles = list(chain.from_iterable(list(reader)))
    Jasp_Angles = [float(x) for x in Jasp_Angles]
    
with open(file_list[1], newline = '') as f:
    reader = csv.reader(f)
    Control_Dist = list(chain.from_iterable(list(reader)))
    Control_Dist = [float(x) for x in Control_Dist]
    
with open(file_list[2], newline = '') as f:
    reader = csv.reader(f)
    Control_Angles = list(chain.from_iterable(list(reader)))
    Control_Angles = [float(x) for x in Control_Angles]
    
with open(file_list[3], newline = '') as f:
    reader = csv.reader(f)
    Control_Length = list(chain.from_iterable(list(reader)))
    Control_Length = [float(x) for x in Control_Length]
    
with open(file_list[4], newline = '') as f:
    reader = csv.reader(f)
    Jasp_Dist = list(chain.from_iterable(list(reader)))
    Jasp_Dist = [float(x) for x in Jasp_Dist]
    
with open(file_list[5], newline = '') as f:
    reader = csv.reader(f)
    Jasp_Length = list(chain.from_iterable(list(reader)))
    Jasp_Length = [float(x) for x in Jasp_Length]
    
All_Dists = Control_Dist + Jasp_Dist
All_Angles = Control_Angles + Jasp_Angles
All_Lengths = Control_Length + Jasp_Length

#Randomize interfilament distances and angles
Control_Angles_Randomized = [random.uniform(0, 90) for x in Control_Angles]
Jasp_Angles_Randomized = [random.uniform(0, 90) for x in Jasp_Angles]
All_Angles_Randomized = [random.uniform(0, 90) for x in All_Angles]

Control_Dist_Randomized = [random.uniform(0, 1000) for x in Control_Dist]
Jasp_Dist_Randomized = [random.uniform(0, 1000) for x in Jasp_Dist]
All_Dists_Randomized = [random.uniform(0, 1000) for x in All_Dists]
"""
#Threshold the interfilament distances
minDist = 0
maxDist = 300

ControlDistFiltered = [x for x in Control_Dist if x > minDist or x < maxDist]
JaspDistFiltered = [x for x in  Jasp_Dist if x > minDist or x < maxDist]
"""

#Plot Control Interfilament distance vs orientation
twod_filename = "Control_InterfilamentDistance_vs_Orientation"
twod_title = "Control Interfilament Distance"

axrange = [[0,300], [0,90]]
fig1, ax1 = plt.subplots()
_, _, _, im = ax1.hist2d(Control_Dist, Control_Angles, bins=50, range=axrange, cmap=my_cmap)
fig1.colorbar(im, ax=ax1, label="4 nm stretches of F-actin per bin")
ax1.set_title(twod_title, wrap=True, pad=10)
ax1.set_yticks(range(0,91,15))
ax1.set_ylabel("Interfilament orientation (º)")
ax1.set_xticks(range(0,301,30))
ax1.set_xlabel("Interfilament distance (nm)")
#fig1.savefig(twod_filename+".png", dpi=600, bbox_inches='tight')

#Plot Jasp interfilament distance vs orientation
twod_filename = "Jasp_InterfilamentDistance_vs_Orientation"
twod_title = "Jasp Interfilament Distance"

axrange = [[0,300], [0,90]]
fig2, ax2 = plt.subplots()
_, _, _, im = ax2.hist2d(Jasp_Dist, Jasp_Angles, bins=50, range=axrange, cmap=my_cmap)
fig2.colorbar(im, ax=ax2, label="4 nm stretches of F-actin per bin")
ax2.set_title(twod_title, wrap=True, pad=10)
ax2.set_yticks(range(0,91,15))
ax2.set_ylabel("Interfilament orientation (º)")
ax2.set_xticks(range(0,301,30))
ax2.set_xlabel("Interfilament distance (nm)")
#fig2.savefig(twod_filename+".png", dpi=600, bbox_inches='tight')

#Plot all interfilament distance vs orientation
twod_filename = "All_InterfilamentDistance_vs_Orientation"
twod_title = "All Interfilament Distance"

axrange = [[0,300], [0,90]]
fig3, ax3 = plt.subplots()
_, _, _, im = ax3.hist2d(All_Dists, All_Angles, bins=50, range=axrange, cmap=my_cmap)
fig3.colorbar(im, ax=ax3, label="4 nm stretches of F-actin per bin")
ax3.set_title(twod_title, wrap=True, pad=10)
ax3.set_yticks(range(0,91,15))
ax3.set_ylabel("Interfilament orientation (º)")
ax3.set_xticks(range(0,301,30))
ax3.set_xlabel("Interfilament distance (nm)")
#fig3.savefig(twod_filename+".png", dpi=600, bbox_inches='tight')

#Plot filament lengths
Figure4 = plt.figure(4, figsize=(4,4))
ax4 = plt.axes([0.3,0.05,0.65,0.9])
ax4.set_ylabel('Filament lengths (nm)',fontsize=30, labelpad=8)
ax4.set_facecolor(cm.Blues(0.05))
ax4.grid(False)
ax4.tick_params('x', top=False, bottom=False, labelsize=0)
ax4.tick_params('y', right=False, labelsize=16)

ax4.set_ylim(-5,max(All_Lengths)+15)
ax4.set_xlim(-0.3,1.3)
ax4_1 = sns.violinplot(data=[Control_Length,Jasp_Length],color=cm.Blues(0.2))
ax4_2 = sns.swarmplot(data=[Control_Length,Jasp_Length],s=3.5,color='b',alpha=0.9)
#Figure4.savefig('BasalActinLengths.png', dpi=600, bbox_inches='tight')

#GET STATISTICS FOR UNTREATED AND JASP FILAMENT LENGTHS
control_Median = np.median(Control_Length)
cont_p75, cont_p25 = np.percentile(Control_Length, [75, 25])

jasp_Median = np.median(Jasp_Length)
jasp_p75, jasp_p25 = np.percentile(Jasp_Length, [75, 25])

filStatistic, filPvalue = st.mannwhitneyu(Control_Length, Jasp_Length, alternative='greater')

#Plot randomized control interfilament distances vs orientations
twod_filename = "Control_InterfilamentDistance_vs_Orientation_Randomized"
twod_title = "Control Interfilament Distance Randomized"

axrange = [[0,300], [0,90]]
fig5, ax5 = plt.subplots()
_, _, _, im = ax5.hist2d(Control_Dist_Randomized, Control_Angles_Randomized, bins=50, range=axrange, cmap=my_cmap)
fig5.colorbar(im, ax=ax5, label="4 nm stretches of F-actin per bin")
ax5.set_title(twod_title, wrap=True, pad=10)
ax5.set_yticks(range(0,91,15))
ax5.set_ylabel("Interfilament orientation (º)")
ax5.set_xticks(range(0,301,30))
ax5.set_xlabel("Interfilament distance (nm)")
#fig5.savefig(twod_filename+".png", dpi=600, bbox_inches='tight')

#Plot randomized jasp interfilament distances vs orientation
twod_filename = "Jasp_InterfilamentDistance_vs_Orientation_Randomized"
twod_title = "Jasp Interfilament Distance Randomized"

axrange = [[0,300], [0,90]]
fig6, ax6 = plt.subplots()
_, _, _, im = ax6.hist2d(Jasp_Dist_Randomized, Jasp_Angles_Randomized, bins=50, range=axrange, cmap=my_cmap)
fig6.colorbar(im, ax=ax6, label="4 nm stretches of F-actin per bin")
ax6.set_title(twod_title, wrap=True, pad=10)
ax6.set_yticks(range(0,91,15))
ax6.set_ylabel("Interfilament orientation (º)")
ax6.set_xticks(range(0,301,30))
ax6.set_xlabel("Interfilament distance (nm)")
#fig6.savefig(twod_filename+".png", dpi=600, bbox_inches='tight')

#Plot randomized all interfilament distances vs orientation
twod_filename = "All_InterfilamentDistance_vs_Orientation_Randomized"
twod_title = "All Interfilament Distance Randomized"

axrange = [[0,300], [0,90]]
fig7, ax7 = plt.subplots()
_, _, _, im = ax7.hist2d(All_Dists_Randomized, All_Angles_Randomized, bins=50, range=axrange, cmap=my_cmap)
fig7.colorbar(im, ax=ax7, label="4 nm stretches of F-actin per bin")
ax7.set_title(twod_title, wrap=True, pad=10)
ax7.set_yticks(range(0,91,15))
ax7.set_ylabel("Interfilament orientation (º)")
ax7.set_xticks(range(0,301,30))
ax7.set_xlabel("Interfilament distance (nm)")
#fig7.savefig(twod_filename+".png", dpi=600, bbox_inches='tight')