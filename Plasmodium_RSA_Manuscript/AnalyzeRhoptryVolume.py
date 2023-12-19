#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 14:16:17 2021

@author: matthewmartinez

Analyze rhoptry volume data. For more detailed info on the code used here, see "AV_measurements2.py"
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

#Hardcode volume data values into 3 lists
AV = [4.37, 3.83, 7.40, 8.46, 7.85, 5.72, 5.14, 6.46, 7.10, 4.36]
Pre = [5.39, 5.45, 7.25, 7.15, 5.56, 4.28, 5.13, 8.39, 4.34]
Fus = [7.14, 4.65, 5.15, 5.93, 6.81, 8.75, 6.34, 5.41, 4.54, 6.73]

Figure1 = plt.figure(1, figsize=(4,4))
ax1 = plt.axes([0.3,0.05,0.8,0.6])
#ax1.set_ylabel('Total rhoptry volume (x1^10 voxels)',fontsize=30, labelpad=8)
#ax1.set_axis_bgcolor(plt.cm.Blues(0.05))
ax1.grid(False)
ax1.tick_params('x', top=False, bottom=False, labelsize=0)
ax1.tick_params('y', left=True, labelsize=16)
ax1.set_ylim(3,9)
ax1.set_xlim(-0.3,0.3)
ax1 = sns.swarmplot(data=[AV, Pre, Fus],s=12,color='b',alpha=0.9)
sns.boxplot(data=[AV, Pre, Fus], color='w', showmeans=True, meanline=True, meanprops={"lw":2,"color":"black"})
