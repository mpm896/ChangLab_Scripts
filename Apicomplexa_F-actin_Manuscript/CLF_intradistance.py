import numpy as np
from matplotlib import pyplot as plt
from sys import argv
import glob
from imodmodel import load_from_modfile
FILTER_DISTANCE = 50
INTERPOLATION_PERIOD = 4

model_folder = "Models/"

#model_folder = argv[1]
file_list = glob.glob("{}/*.mod".format(model_folder))
total_distances = []
total_min_distances = []
for filename in file_list:
	prefix = filename.split("/")[-1].split(".")[0]
	print(prefix)
	model = load_from_modfile(filename)
	#keys = model.get_objects().keys()
	if not "CLF" in filename:
		print("Model {} does not have a CLF object - skipping".format(filename))
		continue
    
	#if len(model.get_object("name").get_contours()) < 2:
        #print("Model {} does not have multiple CLF contours - skipping".format(filename))
		#continue
    
	# if not "microtubule" in keys:
	# 	print("Model {} does not have a microtubule object - skipping".format(filename))
	# 	continue
	# if len(model.get_object("microtubule").get_contours()) == 0:
	# 	print("Model {} does not have any microtubule contours - skipping".format(filename))
	# 	continue
	total_distances_per_model = []
	total_min_distances_per_model = []
	for clf_count, contour in enumerate(model.get_object("name").get_contours()):

		distances = [np.linalg.norm(np.subtract(contour.interpolated_vertices[i+1],contour.interpolated_vertices[i])) for i in range(len(contour.interpolated_vertices)-1)]
		distances_accrued = [sum(distances[0:i]) for i in range(len(distances)+1)]
		# print(len(distances), len(distances_accrued))
		# print(distances_accrued)
		# print(distances)
		# plt.title(septin_count+1)
		# plt.hist(distances, bins=30)
		# plt.show()
		min_distance = 10000
		nearest_clf = None
		distance_vector = [10000,10000,10000]
		# distance_vector=[10000 for i in contour.interpolated_vertices]
		for clf2, clf2_contour in enumerate(model.get_object("name").get_contours()):
			if clf2 == clf_count: # skip self-check
				continue
			vector, distance = contour.contour_nearest_distance(clf2_contour)
			# if distance < min_distance:
			# 	min_distance = distance
			# 	distance_vector = vector
			# 	nearest_septin = septin2
			total_min_distances.append(distance)
			total_min_distances_per_model.append(distance)
			if distance<FILTER_DISTANCE:
				total_distances.extend(vector)
				total_distances_per_model.extend(vector)
				print("CLF Number: ", clf_count+1)
				print("Compared to CLF: ", clf2+1)
				print("Minimum distance: ", distance, "nm")
				# fig, ax = plt.subplots()
				# ax.plot(distances_accrued, distance_vector)
				# ax.set_ylim(0,300)
				# ax.set_ylabel("Distance to nearest septin (nm)")
				# ax.set_xlabel("Distance along septin filament (nm)")
				# ax.set_title("Distance between Septin {} and Septin {}".format(septin_count+1, septin2+1))
				# plt.tight_layout()
				# fig.savefig("septin_{}_{}.svg".format(prefix, septin_count+1))
				# fig2,ax2=plt.subplots()
				# ax2.set_title("Histogram of Distances between Septin {} and Septin {}".format(septin_count+1, septin2+1))
				# ax2.set_ylabel("Count of {}nm stretches of septin".format(INTERPOLATION_PERIOD))
				# ax2.set_xlabel("Septin-MT distance (nm)")
				# ax2.hist(distance_vector,bins=range(0, 100, 4))
				# fig2.savefig("septin_{}_{}_histogram.svg".format(prefix, septin_count+1))
				# print(distance_vector)
		# else:
		# 	print("Septin Number: ", septin_count+1)
		# 	print("No other septins within {}nm".format(FILTER_DISTANCE))
	fig,ax = plt.subplots()
	ax.set_title("Histograms of minimum CLF-CLF distances for all pairs")
	ax.set_ylabel("Count of filament-filament distances".format(INTERPOLATION_PERIOD))
	ax.set_xlabel("CLF-CLF distance (nm)")
	ax.hist(total_min_distances_per_model, bins=range(0, 100, 4))
	#fig.savefig("CLF_total_min_histogram_{}.svg".format(prefix))
	fig3,ax3 = plt.subplots()
	ax3.set_title("Histograms of CLF-CLF distances for all\n CLFs that come within {}nm of another CLF".format(FILTER_DISTANCE))
	ax3.set_ylabel("Count of {}nm stretches of CLF".format(INTERPOLATION_PERIOD))
	ax3.set_xlabel("CLF-CLF distance (nm)")
	ax3.hist(total_distances_per_model, bins=range(0, 100, 4))
	#fig3.savefig("CLF_total_histogram_{}.svg".format(prefix))
	plt.cla()
	plt.clf()

fig,ax = plt.subplots()
ax.set_title("Histograms of minimum CLF-CLF distances for all pairs")
ax.set_ylabel("Count of filament-filament distances".format(INTERPOLATION_PERIOD))
ax.set_xlabel("CLF-CLF distance (nm)")
ax.hist(total_min_distances, bins=range(0, 100, 4))
#fig.savefig("CLF_total_min_histogram.svg")

fig3,ax3 = plt.subplots()
ax3.set_title("Histograms of CLF-CLF distances for all\n septins that come within {}nm of another CLF".format(FILTER_DISTANCE))
ax3.set_ylabel("Count of {}nm stretches of CLF".format(INTERPOLATION_PERIOD))
ax3.set_xlabel("CLF-CLF distance (nm)")
ax3.hist(total_distances, bins=range(0, 100, 2))
#fig3.savefig("CLF_total_histogram.svg")
