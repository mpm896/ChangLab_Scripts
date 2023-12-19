#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:21:49 2022

@author: mattmartinez
Yi-Wei Chang Lab

Script to convert an Amira model into an IMOD model
For example, points along segmented actin filaments in Amira to be converted to
an IMOD model for depiction in ChimeraX

Prior to this (to avoid ridiculous manual parsing) I had to load the Amira XML
file into Microsoft Excel (NOT LibreOffice) and save the "points" and "segments"
into separate files. I couldn't parse using ElementTree.parse because the XML 
file encoding is "System" instead of the normal UTF-8, which caused problems
"""

import glob


pointFiles = glob.glob("*/*Points.csv")
segmentFiles = glob.glob("*/*Segments.csv")

if len(pointFiles) != len(segmentFiles):
    print("Unequal number of Point and Segment files!")
    exit()

for j,file in enumerate(pointFiles):
    key = file.split('/')[0]
    with open(file, 'r') as f:
        points = f.readlines()
        for i in range(len(points)):
            #points[i].decode('utf8').encode('ascii', errors='ignore')
            points[i] = points[i].encode("ascii", "ignore").decode() #Gets rid of hex characters
            points[i] = points[i].strip().split(',')
            points[i] = [float(x) for x in points[i]]
    
    segmentEnds = []
    with open(segmentFiles[j], 'r') as f:
        segmentPoints = f.readlines()
        for i in range(len(segmentPoints)):
            segmentPoints[i] = segmentPoints[i].encode("ascii", "ignore").decode() #Gets rid of hex characters
            segmentPoints[i] = segmentPoints[i].replace('"', '')
            segmentPoints[i] = segmentPoints[i].split(',')
            segmentPoints[i] = [int(x) for x in segmentPoints[i]]
            
            endPoints = [segmentPoints[i][0], segmentPoints[i][-1]]
            segmentEnds.append(endPoints)
    
    contourPoints = []
    for i in range(len(segmentEnds)):
        start = segmentEnds[i][0]
        end = segmentEnds[i][1]
        filamentPoints = points[start:end+1]
        for k in range(len(filamentPoints)):
            filamentPoints[k] = filamentPoints[k][1:4]
        contourPoints.append(filamentPoints)
    
    #Construct the IMOD model file
    fileName = "{}/{}_ActinModel.txt".format(key,key)
    with open(fileName, 'w') as f:
        for i in range(len(contourPoints)):
            for k in range(len(contourPoints[i])):
                line = "{} {} {} {}\n".format(i+1,contourPoints[i][k][0],1440-contourPoints[i][k][1],contourPoints[i][k][2])
                f.write(line)