#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 24 10:00:06 2021

@author: matthewmartinez

Plot data regarding rhoptry morphology class (AV or no AV) vs. membrane integrity vs. treatment
"""
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
import pandas as pd
import seaborn as sns

sns.set()

total = 211
total_intact_AV = 79
total_intact_noAV = 19
total_broken_AV = 65
total_broken_noAV = 48
total_data = [total_intact_AV, total_intact_noAV, total_broken_AV, total_broken_noAV]

total_GlyA= 66
total_GlyA_intact_AV = 26
total_GlyA_intact_noAV = 3
total_GlyA_broken_AV = 21
total_GlyA_broken_noAV = 16
GlyA_data = [total_GlyA_intact_AV, total_GlyA_intact_noAV, total_broken_AV, total_broken_noAV]

total_untreated = 145
total_untreated_intact_AV = 53
total_untreated_intact_noAV = 16
total_untreated_broken_AV = 44
total_untreated_broken_noAV = 32
untreated_data = [total_untreated_intact_AV, total_untreated_intact_noAV, total_untreated_broken_AV, total_untreated_broken_noAV]

df = pd.DataFrame({
        'Treatment': ['GlyA Intact','GlyA Intact','GlyA Intact','GlyA Broken','GlyA Broken','GlyA Broken', 'Untreated Intact','Untreated Intact','Untreated Intact','Untreated Broken','Untreated Broken','Untreated Broken'],
        'AV': ['AV','Prefusion','Fusion','AV','Prefusion','Fusion','AV','Prefusion','Fusion','AV','Prefusion','Fusion'],
        'Value': [26, 1, 2, 21, 7, 9, 53, 6, 10, 44, 18, 14],
        'Percent': [90,3,7,57,19,24,77,9,14,58,23,19]
        })

df2 = pd.DataFrame({
        'Treatment_2': ['GlyA','GlyA','GlyA','Untreated','Untreated','Untreated'],
        'AV_2': ['AV','Prefusion','Fusion','AV','Prefusion','Fusion'],
        'Value_2': [47,8,11,97,24,24],
        'Percent_2': [71,12,17,67,16,16]
        })

df3 = pd.DataFrame({
        'Membrane': ['Intact','Intact','Intact','Broken','Broken','Broken'],
        'AV': ['AV','Prefusion','Fusion','AV','Prefusion','Fusion'],
        'GlyA_Value': [29,1,2,21,7,9],
        'GlyA_Percent': [91,3,6,57,19,24],
        'GlyA_PercentOfTotal': [29,1,2,19,6,8],
        'Untreated_Value': [53,6,10,44,18,14],
        'Untreated_Percent': [77,9,14,58,24,18],
        'Untreated_PercentOfTotal': [52,6,10,39,16,12],
        'Total': [81,7,12,58,22,20]
        })


Figure1 = plt.figure(1,figsize=(4,4))
ax1 = plt.axes([0.3,0.05,0.65,0.9])
ax1.set_ylim(0,100)
ax1 = sns.barplot(x='Treatment_2',y='Percent_2',hue='AV_2',data=df2,ax=ax1)
ax1.legend_.remove()

Figure2 = plt.figure(2, figsize=(10,6))
ax2 = plt.axes([0.3,0.05,0.65,0.9])
ax2 = sns.barplot(x='Treatment',y='Percent',hue='AV',data=df, ax=ax2)
ax2.legend_.remove()

#df3.set_index('Membrane').plot(kind='bar',hue='AV',stacked=True)
Figure3 = plt.figure(3, figsize=(4,4))
ax3 = plt.axes([0.3,0.05,0.65,0.9])
ax3.set_ylim(0,100)
bar1 = sns.barplot(x='Membrane',y='Total',hue='AV',data=df3,hatch='x',ax=ax3)
bar2 = sns.barplot(x='Membrane',y='Untreated_PercentOfTotal',hue='AV',data=df3,ax=ax3)
ax3.legend_.remove()

