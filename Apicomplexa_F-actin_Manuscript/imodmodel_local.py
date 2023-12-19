#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:33:02 2022
Matthew Martinez
Yi-Wei Chang lab

Evaluate distances between filaments of the same kind and different kinds.
Evaluate orientations between filaments of the same kind and different kinds.

Adapted from the Grotjahn Lab "measure models" found here:
    https://github.com/GrotjahnLab/measure_models
"""

import numpy as np
import copy
from scipy.interpolate import splprep, splev
from scipy.stats import binned_statistic
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation as R
from matplotlib import pyplot as plt
import matplotlib
from matplotlib import cm
from sys import argv
import glob
import re

my_cmap = copy.copy(matplotlib.cm.get_cmap('viridis')) # copy the default cmap
my_cmap.set_bad(my_cmap.colors[0])

FILTER_DISTANCE = 10000 #nm, used to find nearby filaments
FILTER_ORIENTATION = 100 
FILTER = None
INTERPOLATE = True
INTERPOLATION_PERIOD = 4 #nm, size chunk of filament to use for analysis
RANDOMIZE = False
objectnames_actin = ["actin", "glideosome actin"]
radii = {"actin": 3.5, "glideosome actin":3.5, "CLF":7}

class IMODModel():
    """Top level model class for data from IMOD"""
    def __init__(self, name, pixel_size, model_range_tuple, offset_tuple):
        self.name               = name
        self.number_of_objects  = 0
        self.pixel_size         = pixel_size
        self.model_range        = model_range_tuple
        self.offset             = offset_tuple
        self.objects            = {}
        
    def add_object(self, imodobject):
        self.number_of_objects  += 1
        self.objects[imodobject.name]   = imodobject
        
    def get_objects(self):
        return self.objects
    
    def get_object(self, name):
        return self.objects[name]
    
    
class IMODObject():
    """Object class for data from IMOD"""
    def __init__(self, name):
        self.name               = name
        self.number_of_contours = 0
        self.contours           = []
        
    def add_contour(self, imodcontour):
        self.number_of_contours += 1
        self.contours.append(imodcontour)
        
    def get_contours(self):
        return self.contours
    
    
class IMODContour():
    """Contour class for data from IMOD based on a list of connected vertices"""
    def __init__(self, vertices, interpolate=True, interpolation_period_nm=INTERPOLATION_PERIOD, randomize=False, tomorange=None):
        vertices = np.array(vertices)
        self.original_vertices  = vertices
        
        vectors = np.diff(vertices, axis=0)
        self.original_vectors = np.append(vectors,vectors[-1:],axis=0)
        assert len(self.original_vectors) == len(self.original_vertices)
        
        self.calculate_length_nm()
        if interpolate:
            self.interpolate_vertices(INTERPOLATION_PERIOD, randomize=randomize, tomorange=tomorange)
            
    def calculate_length_nm(self):
        if len(self.original_vertices) < 2:
            self.length = 0
        else:
            length = 0
            for i in range(len(self.original_vertices)-1):
                length += np.linalg.norm(np.subtract(self.original_vertices[i+1],self.original_vertices[i]))
            self.length = length
        return self.length
    
    def interpolate_vertices(self,interpolation_period_nm, randomize=False, tomorange=None):
        self.interpolation_period_nm = interpolation_period_nm
        if len(self.original_vertices) < 2:
            self.interpolated_vertices = self.original_vertices
        else:
            #TO DO: FIGURE OUT HOW TO MAKE THIS NOT HAVE ROUNDING ERRORS
            x,y,z = zip(*self.original_vertices)
            
            if len(x) > 3:
                tck,u = splprep([x,y,z], k=1)
            else:
                tck,u = splprep([x,y,z], k=1)
                
            u_fine = np.linspace(u[0],u[-1],int(self.length/interpolation_period_nm))
            x_fine,y_fine,z_fine = splev(u_fine,tck)
            
            self.interpolated_vertices =  np.array(
                [x_fine, y_fine, z_fine]).transpose()
            vectors = np.diff(self.interpolated_vertices, axis=0)
            self.interpolated_vectors = np.append(vectors, vectors[-1:], axis=0)
            
            if randomize:
                self.interpolated_vectors,_, self.interpolated_vertices = self.randomize_positions_orientations(self.interpolated_vectors)
            
    def contour_nearest_distance(self, other_contour, interpolated=True):
        if interpolated:
            vertices_self = self.interpolated_vertices
            vertices_other = other_contour.interpolated_vertices
        else:
            vertices_self = self.original_vertices
            vertices_other = other_contour.original_vertices
            
        if len(vertices_self) == 0 or len(vertices_other) == 0:
            return [10000,10000,10000], 10000
        
        distance_vector = [min([np.linalg.norm(np.subtract(i,j)) for j in vertices_other]) for i in vertices_self]
        min_distance = min(distance_vector)
        
        return distance_vector, min_distance
    
    def contour_orientation(self, other_contour, interpolated=True, radius_offset=0):
        """
        This is redundant with def contour_nearest_distance as this also calculates
        distances, but in addition also captures relative orientations
        """

        if interpolated:
            vertices_self = self.interpolated_vertices
            vertices_other = other_contour.interpolated_vertices
            vectors_self = self.interpolated_vectors
            vectors_other = other_contour.interpolated_vectors
        else:
            vertices_self = self.original_vertices
            vertices_other = other_contour.original_vertices
            vectors_self = self.original_vectors
            vectors_other = other_contour.original_vectors
            
        if len(vertices_self) == 0 or len (vertices_other) == 0:
            return [10000, 10000, 10000], 10000, [30,60,90]
        
        dist_matrix = cdist(vertices_self, vertices_other)
        distance_vector = dist_matrix.min(
            axis = 1)-radius_offset
        neighbor_vectors = vectors_other[dist_matrix.argmin(axis=1)]
        # This is a huge performance bottleneck right here - need to vectorize!!!
        relative_orientation = [np.arccos(np.abs(np.dot(vectors_self[i], neighbor_vectors[i]))/(np.linalg.norm(vectors_self[i])*np.linalg.norm(neighbor_vectors[i])))*180/np.pi for i in range(len(vectors_self))]
        min_distance = min(distance_vector)
        
        return distance_vector, min_distance, relative_orientation
    
    def randomize_positions_orientations(self, vector, tomorange = (1000,1000,600), rotrange=[-180, 180], tiltrange=[60,120]):
        """
        Return the same length vector, but with randomized reasonable orientation and fully randomized starting position. 
        Allows constrained rotations to keep general distribution similar.

        Used for generating random distributions with the same overall density. Does not account for excluded volume from organelles, which is standard in PySeg's similar routines.

        Randomizes rot and tilt - because of the way euler angles are selected, large tilt ranges may lead to unacceptably nonuniform sampling.

        Currently is only applied to interpolated vectors, as later interpolation is messed up by wrapping.
        """
        ### Calculate and correct current orientations ###
        #tomorange = [i/10 for i in tomorange]
        initial_vector = vector[0]
        initial_theta = np.arctan2(initial_vector[1],initial_vector[0])*180/np.pi
        initial_phi = np.arccos(initial_vector[2]/np.linalg.norm(initial_vector))*180/np.pi
        correction_rotation = R.from_euler("ZYZ", [initial_theta, initial_phi,200], degrees=True)
        correction_rotation = correction_rotation.inv()
        corrected_vector = correction_rotation.apply(vector)
        
        ### Randomize angles and positions ###
        initial = np.random.rand(3)*tomorange # rand is from 0 to 1
        rot = np.random.rand()*(rotrange[1]-rotrange[0])+rotrange[0]
        tilt = np.random.rand()*(tiltrange[1]-tiltrange[0])+tiltrange[0]
        rotation = R.from_euler("ZY", [rot, tilt], degrees=True)
        rotated_vector = rotation.apply(corrected_vector,)
        rotated_positions = np.cumsum(rotated_vector, axis=0)+initial
        rotated_positions_wrapped = rotated_positions
        
        for i in range(3):
            rotated_positions_wrapped[:][i] = rotated_positions[:][i] % tomorange[i]
            
        return rotated_vector, rotated_positions, rotated_positions_wrapped
        


def load_from_modfile(filename):
    import subprocess
    
    #Initialize model variables
    asciifile           = subprocess.run(["imodinfo", "-a", filename], check=True, capture_output=True).stdout.decode()
    model               = None
    pixel_size          = 1
    model_range         = None
    objects             = None
    contours            = None 
    number_of_objects   = None
    number_of_contours  = None
    offsets             = None 
    angles              = None
    
    lines = asciifile.split("\n")
    iterator = iter(lines)
    oldvertex = (0,0,0)
    
    #Iterate through the lines and initialize model/objects/contours
    for line in iterator:
        if line.startswith("#"):
            continue
        if line == '':
            continue
        if line.startswith("imod"):
            number_of_objects   = int(line.split()[-1])
        if line.startswith("max"):
            model_range         = (int(i) for i in line.split()[1:4])
        if line.startswith("offsets"):
            offsets             = (float(i) for i in line.split()[1:4])
        if line.startswith("angles"):
            angles              = (float(i) for i in line.split()[1:4])
        if line.startswith("pixsize"):
            pixel_size          = float(line.split()[-1])
        if line.startswith("units"):
            #INITIALIZE THE MODEL
            model = IMODModel(filename, pixel_size, model_range, offsets)
        if line.startswith("object"):
            object_number, number_of_contours, number_of_meshes = (int(i) for i in line.split()[1:4])
            name = next(iterator).split(' ', 1)[1] #removes "name" and saves the entire object name
            if "knee" in name:
                continue
            if name == '':
                name = str(object_number) #Names object by number if there are no names
            
            #INITIALIZE EACH OBJECT IN THE MODEL
            newobject = IMODObject(name)
            model.add_object(newobject)
            
            newline = next(iterator)
            while newline != "":
                if newline.startswith("contour"):
                    contour_number, surface, n_vertices = (int(i) for i in newline.split()[1:4])
                    
                    vertices = []
                    for i in range(n_vertices):
                        newline = next(iterator)
                        vertex = tuple([float(i)*pixel_size for i in newline.split()])
                        if not oldvertex == vertex:
                            vertices.append(vertex)
                            oldvertex = vertex
                    
                    #INITIALIZE EACH CONTOUR IN THE OBJECT
                    newobject.add_contour(IMODContour(vertices, interpolate=True, interpolation_period_nm=INTERPOLATION_PERIOD, randomize=RANDOMIZE))
                newline = next(iterator)
                
    return model


def Fil_Fil_Distances():
    """
    Gets interfilament distances between CLFs and all F-actin filaments,
    regardless of F-actin location (conoid vs. IMC-PM space)
    """
    import re
    
    #model_folder           = argv[1]
    model_folder            = "Models/"
    file_list               = glob.glob("{}/*.mod".format(model_folder))
    tomo_ID_list            = (set([file_list[i].split("/")[-1].split("_")[0] for i in range(len(file_list))])) #Gets list of tomogram models, assuming different filament types were made in different model files
    print(tomo_ID_list)
    
    total_distances         = []
    total_distances_accrued = []
    sum_distance            = 0
    sum_distance_nearby     = 0
    near_actin_count        = 0 
    actin_total_count       = 0
    clf_total_count         = 0
    
    for tomo_ID in tomo_ID_list:
        mod_files = glob.glob("{}/{}*.mod".format(model_folder, tomo_ID))
        print(mod_files)
        
        near_actin_count = 0
        near_actin_count_file = 0
        sum_distance_per_file = 0
        sum_distance_nearby_per_file = 0
    
        print(tomo_ID)
        
        #Determine which actin models are present in each tomogram
        #INITIALIZE CLF AND F-ACTIN MODELS IF MODEL MADE FOR TOMOGRAM
        model_actin = None
        model_clf = None
        model_actin = [load_from_modfile(mod_files[i]) for i in range(len(mod_files)) if re.search('actin', mod_files[i], flags=re.IGNORECASE)]
        model_clf =  [load_from_modfile(mod_files[i]) for i in range(len(mod_files)) if re.search('clf', mod_files[i], flags=re.IGNORECASE)]
        
        #Run this code of the model of F-actin exists for this tomogram
        if model_actin: 
            model_actin = model_actin[0]
            keys_actin = model_actin.get_objects().keys()
            actin_keys = ["actin", "glideosome actin"]
            match = [c for c in actin_keys if c in list(keys_actin)]
            print(list(keys_actin))
    
            if len(match) == 0:
                print("No actin in this model")
                continue
                    
            actin_count_per_file = 0
            glideosome_actin_per_file = 0
        
            for key in match:
                if key == "actin":
                    actin_count_per_file = len(model_actin.get_object("actin").get_contours())
                if key == "glideosome actin":
                    glideosome_actin_per_file = len(model_actin.get_object("glideosome actin").get_contours())
                    actin_total_count += actin_count_per_file + glideosome_actin_per_file
            print("Actin count: " + str(actin_count_per_file) + " " + str(glideosome_actin_per_file))
        
        #Run this code of the model of the CLFs exists for this tomogram
        if model_clf:
            model_clf = model_clf[0]
            keys_clf = model_clf.get_objects().keys()
             
            if len(model_clf.get_objects()) == 0:
                 print("No CLFs in this model")
                 continue
         
            clf_count_per_object = 0
            for key in list(keys_clf):
                 clf_count_per_object += len(model_clf.get_object(key).get_contours())
                 clf_count_per_file = clf_count_per_object
                 clf_total_count += clf_count_per_file
            print("CLF count:" + str(clf_count_per_file))  
            
            #CALCULATE SPATIAL RELATIONSHIP BETWEEN CLFs AND F-ACTIN
            #Need to iterate separately for models with multiple CLF objects
            #(one contour per object) amd for models with one object that contain
            #all the CLF contours, because I messed up these models
            total_distance_per_model = []
            if len(model_clf.get_objects()) == 1:
                for clf_count, contour in enumerate(model_clf.get_object(list(keys_clf)[0]).get_contours()):
                    distances = [np.linalg.norm(np.subtract(contour.interpolated_vertices[i+1],contour.interpolated_vertices[i])) for i in range(len(contour.interpolated_vertices)-1)]
                    distances_accrued = [sum(distances[0:i]) for i in range(len(distances)+1)]
                    sum_distance_per_clf = distances_accrued[-1]
                    sum_distance_per_file += sum_distance_per_clf
                    min_distance = 10000
                    nearest_actin = None
                    distance_vector = [10000,10000,10000]
                    
                    if model_actin:
                        for j in range(len(match)):
                            for actin_position, actin_contour in enumerate(model_actin.get_object(match[j]).get_contours()):
                                vector, distance = contour.contour_nearest_distance(actin_contour)
                                if distance < min_distance:
                                    min_distance = distance
                                    distance_vector = vector
                                    nearest_actin = actin_position
                                    
                    if min_distance < FILTER_DISTANCE:
                        near_actin_count += 1
                        near_actin_count_file += 1
                        sum_distance_nearby_per_clf = sum([k for j,k in enumerate(distances) if distance_vector[j]<FILTER_DISTANCE])
                        sum_distance_nearby_per_file += sum_distance_per_clf
                        total_distances.extend(distance_vector)
                        total_distance_per_model.extend(distance_vector)
                        total_distances_accrued.extend(distances_accrued)
                        
                        #Print basic analysis values
                        """
                        print("CLF number: ", clf_count+1)
                        print("Nearest F-actin: ", nearest_actin+1)
                        print("Minimum distance: ", min_distance, "nm")
                        print("Total length of CLF: ", sum_distance_per_clf, "nm")
                        print("Length of CLF within {}nm of F-actin {}: {}nm".format(FILTER_DISTANCE, nearest_actin+1, sum_distance_nearby_per_clf))
                        """
                        
                        #Plot results
                        fig, ax = plt.subplots()
                        ax.plot(distances_accrued, distance_vector)
                        ax.set_ylim(0,300)
                        ax.set_ylabel("Distance to an F-actin (nm)")
                        ax.set_xlabel("Distance along CLF (nm)")
                        ax.set_title("Distance between CLF {} and F-actin {}".format(clf_count+1, nearest_actin+1))
                        plt.tight_layout()
                        
                        
                        fig2, ax2 = plt.subplots()
                        ax2.hist(distance_vector,bins=range(0, 100, 4))
                        ax2.set_xticks([0,10,20,30,40,50,60,70,80,90,100])
                        ax2.set_xlabel("Distance (nm)")
                        ax2.set_title("Histogram of Distances between CLF {} and F-actin {}".format(clf_count+1, nearest_actin+1))
                        ax2.set_ylabel("Count of {}nm stretches of CLF".format(INTERPOLATION_PERIOD))
                        ax2.set_xlabel("CLF-Actin distance (nm)")
                        plt.clf()
                    else:
                        sum_distance_nearby_per_clf = 0
                        """
                        print("CLF number: ", clf_count+1)
                        print("No F-actin within {}nm".format(FILTER_DISTANCE))
                        """
                        
                fig3, ax3 = plt.subplots()
                ax3.set_title("Histograms of CLF-Actin distances for all\n CLFs that come within {}nm of an F-actin in tomo {}".format(FILTER_DISTANCE, tomo_ID))
                ax3.set_ylabel("Count of {}nm stretches of CLF".format(INTERPOLATION_PERIOD))
                ax3.set_xlabel("CLF-Actin distance (nm)")
                ax3.hist(total_distance_per_model, bins=range(0, 100, 4))
                ax3.set_xticks([0,10,20,30,40,50,60,70,80,90,100])
                ax3.set_xlabel("Distance (nm)")
                #fig3.savefig("Total_histogram_{}.svg".format(prefix))
                #plt.close('all')
                sum_distance += sum_distance_per_file
                sum_distance_nearby += sum_distance_nearby_per_file
                        
            if len(model_clf.get_objects()) > 1:
                clf_count = 0
                for j in range(len(model_clf.get_objects())):
                    contour = (model_clf.get_object(list(keys_clf)[j]).get_contours())[0]
                    distances = [np.linalg.norm(np.subtract(contour.interpolated_vertices[k+1],contour.interpolated_vertices[k])) for k in range(len(contour.interpolated_vertices)-1)]
                    distances_accrued = [sum(distances[0:k]) for k in range(len(distances)+1)]
                    sum_distance_per_clf = distances_accrued[-1]
                    sum_distance_per_file += sum_distance_per_clf
                    min_distance = 10000
                    nearest_actin = None
                    distance_vector = [10000,10000,10000]
                    
                    if model_actin:
                        for k in range(len(match)):
                            for actin_position, actin_contour in enumerate(model_actin.get_object(match[k]).get_contours()):
                                vector, distance = contour.contour_nearest_distance(actin_contour)
                                if distance < min_distance:
                                    min_distance = distance
                                    distance_vector = vector
                                    nearest_actin = actin_position
                                    
                    if min_distance < FILTER_DISTANCE:
                        near_actin_count += 1
                        near_actin_count_file += 1
                        sum_distance_nearby_per_clf = sum([k for l,k in enumerate(distances) if distance_vector[l]<FILTER_DISTANCE])
                        sum_distance_nearby_per_file += sum_distance_per_clf
                        total_distances.extend(distance_vector)
                        total_distance_per_model.extend(distance_vector)
                        total_distances_accrued.extend(distances_accrued)
                        
                        #Print basic analysis values
                        """
                        print("CLF number: ", clf_count+1)
                        print("Nearest F-actin: ", nearest_actin+1)
                        print("Minimum distance: ", min_distance, "nm")
                        print("Total length of CLF: ", sum_distance_per_clf, "nm")
                        print("Length of CLF within {}nm of F-actin {}: {}nm".format(FILTER_DISTANCE, nearest_actin+1, sum_distance_nearby_per_clf))
                        """
                        
                        #Plot results
                        fig, ax = plt.subplots()
                        ax.plot(distances_accrued, distance_vector)
                        ax.set_ylim(0,300)
                        ax.set_ylabel("Distance to an F-actin (nm)")
                        ax.set_xlabel("Distance along CLF (nm)")
                        ax.set_title("Distance between CLF {} and F-actin {}".format(clf_count+1, nearest_actin+1))
                        plt.tight_layout()
                        
                        
                        fig2, ax2 = plt.subplots()
                        ax2.hist(distance_vector,bins=range(0, 100, 4))
                        ax2.set_xticks([0,10,20,30,40,50,60,70,80,90,100])
                        ax2.set_xlabel("Distance (nm)")
                        ax2.set_title("Histogram of Distances between CLF {} and F-actin {}".format(clf_count+1, nearest_actin+1))
                        ax2.set_ylabel("Count of {}nm stretches of CLF".format(INTERPOLATION_PERIOD))
                        ax2.set_xlabel("CLF-Actin distance (nm)")
                        plt.clf()
                    else:
                        sum_distance_nearby_per_clf = 0
                        """
                        print("CLF number: ", clf_count+1)
                        print("No F-actin within {}nm".format(FILTER_DISTANCE))
                        """
                    clf_count += 1
                
                
                fig3, ax3 = plt.subplots()
                ax3.set_title("Histograms of CLF-Actin distances for all\n CLFs that come within {}nm of an F-actin in tomo {}".format(FILTER_DISTANCE, tomo_ID))
                ax3.set_ylabel("Count of {}nm stretches of CLF".format(INTERPOLATION_PERIOD))
                ax3.set_xlabel("CLF-Actin distance (nm)")
                ax3.hist(total_distance_per_model, bins=range(0, 100, 4))
                ax3.set_xticks([0,10,20,30,40,50,60,70,80,90,100])
                ax3.set_xlabel("Distance (nm)")
                #fig3.savefig("Total_histogram_{}.svg".format(prefix))
                #plt.close('all')
                    
                sum_distance += sum_distance_per_file
                sum_distance_nearby += sum_distance_nearby_per_file
              
    fig3, ax3 = plt.subplots()
    ax3.set_title("Histograms of CLF-Actin distances for all\n CLFs that come within {}nm of an F-actin".format(FILTER_DISTANCE*2))
    ax3.set_ylabel("Count of {}nm stretches of CLF".format(INTERPOLATION_PERIOD))
    ax3.set_xlabel("CLF-Actin distance (nm)")
    ax3.hist(total_distances, bins=range(0, FILTER_DISTANCE*2, 4))
    ax3.set_xticks(0, FILTER_DISTANCE*2, 10)
    ax3.set_xlabel("Distance (nm)")
    #fig3.savefig("Total_histogram_{}.svg".format(prefix))
    
    binned_distances = binned_statistic(total_distances_accrued, total_distances, statistic='mean', bins=20, range=(0,200))
    
    fig4, ax4 = plt.subplots()
    ax4.plot(binned_distances.bin_edges[0:-1], binned_distances.statistic)
    ax4.set_ylim(0,300)
    ax4.set_xlim(0,200)
    ax4.set_ylabel("Distance to an F-actin (nm)")
    ax4.set_xlabel("Distance along CLF (nm)")
    ax4.set_title("Distance between all CLFs and F-actins")
    plt.tight_layout()
                
"""
WITH WHERE I'M LEAVING OFF, NEED TO SET CONDITIONALS TO ANALYZE NON-PCR ASSOCIATED
F-ACTIN WITH RESPECT TO CLFs, OR PROVIDE A CUTOFF TO ONLY ANALYZE F-ACTIN 
SEGMENTS BELOW THE APR
"""                   


def distances_orientations(imodmodel, key_1, imodmodel2=None, key_2=None, randomize=False, radius_offset=0, filter=FILTER, individual_plots=False):
    """
    Calculates interfilament distances and orientations. Can be within one type
    of filament (i.e. just between the CLFs) or between two different filament 
    types (i.e. CLFs and F-actin) depending on the models provides
    """
    full_actin_comparison = False
    if key_1 == "CLF" and (not key_2):
       key_2 = key_1
       self_comparison = True
    elif not key_2:
        key_2 = key_1
        self_comparison = True
    elif (key_1 == "actin" or key_2 == "actin") and (key_1 == "glideosome actin" or key_2 == "glideosome actin"):
        full_actin_comparison = True
        self_comparison = True
    else:
        self_comparison = False
    
    min_distance    = 100000
    distances       = []
    orientations    = []
    
    if full_actin_comparison:
        actin_contours = imodmodel.get_object(key_1).get_contours() + imodmodel.get_object(key_2).get_contours()
        for count, contour in enumerate(actin_contours):
            min_distance        = 100000
            nearest_contour     = None
            orientation_vector  = []
            distance_vector     = []
            
            for count_2, contour_2 in enumerate(imodmodel.get_object(key_2).get_contours()):
                if self_comparison and count_2 == count:
                    continue
                vector, distance, orientation = contour.contour_orientation(
                    contour_2, radius_offset=radius_offset)
                
                if distance < min_distance:
                    min_distance = distance
                    distance_vector = vector
                    nearest_contour = count_2
                    orientation_vector = orientation
            
            if filter and min_distance > filter:
                continue
            distances.append(distance_vector)
            orientations.append(orientation_vector)
        distances = np.concatenate(distances, axis=None)
        orientations = np.concatenate(orientations, axis=None)
           
        distances_post = []
        orientations_post = []
        for i in range(len(distances)):
            if distances[i] > 0:
                distances_post.append(distances[i])
                orientations_post.append(orientations[i])
        distances_post = np.concatenate(distances_post, axis=None)
        orientations_post = np.concatenate(orientations_post, axis=None) 
    else:    
        for count, contour in enumerate(imodmodel.get_object(key_1).get_contours()):
            min_distance        = 100000
            nearest_contour     = None
            orientation_vector  = []
            distance_vector     = []
            
            
            #Do this block if there are 2 different models (i.e. comparing actin to CLF)
            if imodmodel2:
                for count_2, contour_2 in enumerate(imodmodel2.get_object(key_2).get_contours()):
                    if self_comparison and count_2 == count:
                        continue
                    vector, distance, orientation = contour.contour_orientation(
                        contour_2, radius_offset=radius_offset)
                    
                    if distance < min_distance:
                        min_distance = distance
                        distance_vector = vector
                        nearest_contour = count_2
                        orientation_vector = orientation
                                
                if filter and min_distance > filter:
                    continue
                distances.append(distance_vector)
                orientations.append(orientation_vector)
            #Do this block if it is a self comparison among CLFs or one type of F-actin
            elif (key_1 == "CLF") or (key_1 != "CLF" and self_comparison):
                for count_2, contour_2 in enumerate(imodmodel.get_object(key_2).get_contours()):
                    if self_comparison and count_2 == count:
                        continue
                    vector, distance, orientation = contour.contour_orientation(
                        contour_2, radius_offset=radius_offset)
                    
                    if distance < min_distance:
                        min_distance = distance
                        distance_vector = vector
                        nearest_contour = count_2
                        orientation_vector = orientation
                                
                if filter and min_distance > filter:
                    continue
                distances.append(distance_vector)
                orientations.append(orientation_vector)
            #Do this block if self comparison of different types of F-actin       
            else:
                for count_2, contour_2 in enumerate(imodmodel.get_object(key_2).get_contours()):
                    if self_comparison and count_2 == count:
                        continue
                    vector, distance, orientation = contour.contour_orientation(
                        contour_2, radius_offset=radius_offset)
                    
                    if distance < min_distance:
                        min_distance = distance
                        distance_vector = vector
                        nearest_contour = count_2
                        orientation_vector = orientation
                
                if filter and min_distance > filter:
                    continue
                distances.append(distance_vector)
                orientations.append(orientation_vector)

        
        if individual_plots and nearest_contour:
            makeplots(imodmodel.name, distances, orientations, key_1+"_"+count, key_2+"_"+nearest_contour, randomize, filter=None, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)
    
        distances = np.concatenate(distances, axis=None)
        orientations = np.concatenate(orientations, axis=None)
        
        distances_post = []
        orientations_post = []
        for i in range(len(distances)):
            if distances[i] > 0:
                distances_post.append(distances[i])
                orientations_post.append(orientations[i])
        distances_post = np.concatenate(distances_post, axis=None)
        orientations_post = np.concatenate(orientations_post, axis=None)
    
    return distances_post, orientations_post
    
def makeplots(tomo_ID, distances, orientations, key1, key2, randomize, filter, interpolate, interpolation_period, format=".svg"):
    """
    Plots data on filament orientations
    """
    #Distances histogram
    distances_filename      = f"{tomo_ID}_{key1}_{key2}_distances"
    distances_title         = f"{tomo_ID} {key1}-{key2} distance histogram"
    if filter:
        distances_filename  += (f"_{filter}nm")
        distances_title     += (f" for {key1} filaments within {filter}nm of any {key2}")
    if randomize:
        distances_filename  += ("_randomized")
        distances_title     += (" (randomized positions and orientations)")
        
    distances_yaxis         = f"Count of {interpolation_period}nm stretches of {key1}"
    distances_xaxis         = f"{key1}-{key2} distance (nm)"
    fig, ax = plt.subplots()
    ax.set_title(distances_title, wrap=True, pad=10)
    ax.set_ylabel(distances_yaxis)
    ax.set_xlabel(distances_xaxis)
    ax.hist(distances, bins=range(0, 201, 4))
    ax.set_xlim(0,200)
    ax.set_xticks(range(0,201,20))
    #fig.savefig(distances_filename+format)
    
    
    #Orientations
    orientations_filename   = f"{tomo_ID}_{key1}_{key2}_orientations"
    orientations_title      = f"{tomo_ID} {key1}-{key2} orientation histogram"
    if filter:
        orientations_filename += f"_{filter}nm"
        orientations_title  += f" for {key1} filaments within {filter}nm of any {key2}"
    if randomize:
        orientations_filename += "_randomized"
        orientations_title  += " (randomized positions and orientations)"
    
    orientations_yaxis      = f"Count of {interpolation_period}nm stretches of {key1}"
    orientations_xaxis      = f"{key1}-{key2} orientation (º)"  
    fig2, ax2 = plt.subplots()
    ax2.set_title(orientations_title, wrap=True, pad=10)
    ax2.set_ylabel(orientations_yaxis)
    ax2.set_xlabel(orientations_xaxis)
    ax2.hist(distances, bins=range(0, 91, 5))
    ax2.set_xticks(range(0,91,15))
    ax2.set_xlim(0,90)
    #fig2.savefig(orientations_filename+format)
    
    
    #Distances vs Orientations
    twod_filename           = f"{tomo_ID}_{key1}_{key2}_2d_hist"
    twod_title              = f"{tomo_ID} {key1}-{key2} orientation vs distance heatmap"
    if filter:
        twod_filename       += f"_{filter}nm"
        twod_title          += f" for {key1} filaments within {filter}nm of any {key2}"
    if randomize:
        twod_filename       += ("_randomized")
        twod_title          += (" (randomized positions and orientations)")
        
    axrange = [[0,200], [0,90]]
    fig3, ax3 = plt.subplots()
    _, _, _, im = ax3.hist2d(distances, orientations, bins=30, range=axrange, cmap=my_cmap)
    fig3.colorbar(im, ax=ax3, label=f"{interpolation_period}nm stretches of {key1} per bin")
    ax3.set_title(twod_title, wrap=True, pad=10)
    ax3.set_yticks(range(0,91,15))
    ax3.set_ylabel(f"{key1}-{key2} orientation (º)")
    ax3.set_xticks(range(0,201,20))
    ax3.set_xlabel(f"{key1}-{key2} distance (nm)")
    #fig3.savefig(twod_filename+format)
    
    fig4, ax4 = plt.subplots()
    _, _, _, im = ax4.hist2d(distances, orientations, bins=30, range=axrange, norm=matplotlib.colors.LogNorm(), cmap=my_cmap)
    fig4.colorbar(im, ax=ax4, label=f"{interpolation_period}nm stretches of {key1} per bin")
    ax4.set_title(twod_title, wrap=True, pad=10)
    ax4.set_yticks(range(0,91,15))
    ax4.set_ylabel(f"{key1}-{key2} orientation (º)")
    ax4.set_xticks(range(0,201,20))
    ax4.set_xlabel(f"{key1}-{key2} distance (nm)")
    fig4.savefig(twod_filename+"_log"+format)
    
    
actin_actin_distances = []
actin_actin_orientations = []

PCR_PCR_distances = []
PCR_PCR_orientations = []

clf_actin_distances = []
clf_actin_orientations = []

clf_PCR_distances = []
clf_PCR_orientations = []

clf_glid_distances = []
clf_glid_orientations = []

clf_clf_distances = []
clf_clf_orientations = []


#model_folder           = argv[1]
model_folder            = "Models/"
file_list               = glob.glob("{}/*.mod".format(model_folder))
tomo_ID_list            = (set([file_list[i].split("/")[-1].split("_")[0] for i in range(len(file_list))])) #Gets list of tomogram models, assuming different filament types were made in different model files
print(tomo_ID_list)

for tomo_ID in tomo_ID_list:
    mod_files = glob.glob("{}/{}*.mod".format(model_folder, tomo_ID))
    print(mod_files)
    
    near_actin_count = 0
    near_actin_count_file = 0
    sum_distance_per_file = 0
    sum_distance_nearby_per_file = 0

    print(tomo_ID)
    model_actin = None
    model_clf = None
    
    #Don't proceed if the tomogram doesn't have both an actin and CLF model
    if len(mod_files) < 2:
        for file in mod_files:
            if not 'actin' in file.lower():
                print("{} doesn't have an associated actin model".format(tomo_ID))
            if not 'clf' in file.lower():
                print("{} doesn't have an associated CLF model".format(tomo_ID))
        continue
    
    #Determine which actin models are present in each tomogram
    #INITIALIZE CLF AND F-ACTIN MODELS IF MODEL MADE FOR TOMOGRAM
    model_actin_tomo = [load_from_modfile(mod_files[i]) for i in range(len(mod_files)) if re.search('actin', mod_files[i], flags=re.IGNORECASE)]
    print("Loaded actin model")
    if len(model_actin_tomo) > 0:
        model_actin = model_actin_tomo[0]
    model_clf_tomo =  [load_from_modfile(mod_files[i]) for i in range(len(mod_files)) if re.search('clf', mod_files[i], flags=re.IGNORECASE)]
    print("Loaded CLF model")
    if len(model_clf_tomo) > 0:
        model_clf = model_clf_tomo[0]
    

    
    if model_clf:
        keys_model_clf = model_clf.get_objects().keys()
        key_clf = list(keys_model_clf)[0]
    if model_actin:
        keys_model_actin = model_actin.get_objects().keys()
        keys_actin = [c for c in objectnames_actin if c in list(keys_model_actin)] #Match only the desired actin objects, "actin" and "glideosome actin"
    
    if len(keys_actin) == 1:
        key = keys_actin[0]
        
        #Get spatial info on just the actin component (models with just 1 actin component)
        radius = radii["actin"]*2
        distance_actin, orientation_actin = distances_orientations(model_actin, key, randomize=RANDOMIZE, radius_offset=radius, filter=FILTER)
        actin_actin_distances.append(distance_actin)
        actin_actin_orientations.append(orientation_actin)
        
        
        #Get spatial info on just the CLF components
        radius = radii["CLF"]*2
        distance_clf, orientation_clf = distances_orientations(model_clf, key_clf, randomize=RANDOMIZE, radius_offset=radius, filter=FILTER)
        clf_clf_distances.append(distance_clf)
        clf_clf_orientations.append(orientation_clf)
        
        #Get spatial info on F-actin vs. CLFs
        radius = radii["CLF"]+radii["actin"]
        distance, orientation = distances_orientations(model_clf, key_clf, imodmodel2=model_actin, key_2=key, randomize=RANDOMIZE, radius_offset=radius, filter=FILTER)
        clf_actin_distances.append(distance)
        clf_actin_orientations.append(orientation)
        
        if len(distance) > 0 and min(distance) < 200:
            makeplots(tomo_ID, distance_actin, orientation_actin, key, key, randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)
            makeplots(tomo_ID, distance_clf, orientation_clf, key_clf, key_clf, randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)
            makeplots(tomo_ID, distance, orientation, key_clf, key, randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)
            
    if len(keys_actin) == 2:
        
        #Get spatial info on just the PCR-actin component
        radius = radii["actin"]*2
        distance_PCR_actin, orientation_PCR_actin = distances_orientations(model_actin, "actin", randomize=RANDOMIZE, radius_offset=radius, filter=FILTER)
        PCR_PCR_distances.append(distance_PCR_actin)
        PCR_PCR_orientations.append(orientation_PCR_actin)
        
        #Get spatial info on the actin component
        radius = radii["actin"]*2
        distance_actin, orientation_actin = distances_orientations(model_actin, "actin", key_2="glideosome actin", randomize=RANDOMIZE, radius_offset=radius, filter=FILTER)
        actin_actin_distances.append(distance_actin)
        actin_actin_orientations.append(orientation_actin)
        
        #Get spatial info on just the CLF components
        radius = radii["CLF"]*2
        distance_clf, orientation_clf = distances_orientations(model_clf, "CLF", randomize=RANDOMIZE, radius_offset=radius, filter=FILTER)
        clf_clf_distances.append(distance_clf)
        clf_clf_orientations.append(orientation_clf)
        
        #Get spatial info on CLF vs PCR actin
        radius = radii["CLF"]+radii["actin"]
        distance_clf_PCR, orientation_clf_PCR = distances_orientations(model_clf, "CLF", imodmodel2=model_actin, key_2="actin", randomize=RANDOMIZE, radius_offset=radius, filter=FILTER)
        clf_PCR_distances.append(distance_clf_PCR)
        clf_PCR_orientations.append(orientation_clf_PCR)
        clf_actin_distances.append(distance_clf_PCR)
        clf_actin_orientations.append(orientation_clf_PCR)
        
        #Get spatial info on CLF vs glideosome actin
        radius = radii["CLF"]+radii["actin"]
        distance_clf_glid, orientation_clf_glid = distances_orientations(model_clf, "CLF", imodmodel2=model_actin, key_2="glideosome actin", randomize=RANDOMIZE, radius_offset=radius, filter=FILTER)
        clf_glid_distances.append(distance_clf_glid)
        clf_glid_orientations.append(orientation_clf_glid)
        clf_actin_distances.append(distance_clf_glid)
        clf_actin_orientations.append(orientation_clf_glid)
        
        if len(distance_clf) > 0 and min(distance_clf) < 200:
            #makeplots(tomo_ID, distance_actin, orientation_actin, "actin", "actin", randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)
            makeplots(tomo_ID, distance_clf, orientation_clf, "CLF", "CLF", randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)
            #makeplots(tomo_ID, distance_clf_PCR, orientation_clf_PCR, "CLF", "PCR actin", randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)
            #makeplots(tomo_ID, distance_clf_glid, orientation_clf_glid, "CLF", "Glideosome actin", randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)
"""        
actin_actin_distances = np.concatenate(actin_actin_distances, axis=None)
actin_actin_orientations = np.concatenate(actin_actin_orientations, axis=None)
makeplots("Total", actin_actin_distances, actin_actin_orientations, "F-actin", "F-actin", randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)

PCR_PCR_distances = np.concatenate(PCR_PCR_distances)
PCR_PCR_orientations = np.concatenate(PCR_PCR_orientations)
makeplots("Total", PCR_PCR_distances, PCR_PCR_orientations, "PCR F-actin", "PCR F-actin", randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)

clf_actin_distances = np.concatenate(clf_actin_distances)
clf_actin_orientations = np.concatenate(clf_actin_orientations)
makeplots("Total", clf_actin_distances, clf_actin_orientations, "CLF", "F-actin", randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)

clf_PCR_distances = np.concatenate(clf_PCR_distances)
clf_PCR_orientations = np.concatenate(clf_PCR_orientations)
makeplots("Total", clf_PCR_distances, clf_PCR_orientations, "CLF", "PCR F-actin", randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)

clf_glid_distances = np.concatenate(clf_glid_distances)
clf_glid_orientations = np.concatenate(clf_glid_orientations)
makeplots("Total", clf_glid_distances, clf_glid_orientations, "CLF", "Glideosome actin", randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)
"""
clf_clf_distances = np.concatenate(clf_clf_distances)
clf_clf_orientations = np.concatenate(clf_clf_orientations)
makeplots("Total", clf_clf_distances, clf_clf_orientations, "CLF", "CLF", randomize=RANDOMIZE, filter=FILTER, interpolate=INTERPOLATE, interpolation_period=INTERPOLATION_PERIOD)
    
    
    
    
"""
WHERE I'M LEAVING OFF, FOLLOW THE CODE AT THE BOTTOM OF "amira_filament_comparator.py"
TO PERFORM ALL ANALYSIS AND PLOTTING
"""
    
    
    
    
    
    
