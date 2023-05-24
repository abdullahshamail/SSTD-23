# Import all libraries

import numpy as np
from sklearn.cluster import DBSCAN
import pandas as pd
import matplotlib.pyplot as plt
from utilities import hydrogenBondCheckGivenIndices2, calculateAtomTypes, printConvoys, consecutive_true_indices, combine_ranges
from CMC import CMC
from hyperparameters import Hyperparameters
import time


#define the hyperparameters
hyperParameters = Hyperparameters(m = 50, k = 5, epsilon = 3.5, t1 = 1, t2 = 50, threshold=4)
dataSize = 500

if __name__ == '__main__':

    fileLocation = r"...\molecularInformation.csv" #location for information of atoms

    filePath = r"...\xyzCoords.npy"

    data = np.load(filePath)
    data = data[:dataSize]
    data1 = list()
    for x in range(data.shape[1]):
        data1.append(data[:,x,0:3].tolist())


    # Clustering using DBSCAN

    totalTime = 0

    t1 = time.time()

    clustering_clf = DBSCAN(eps=hyperParameters.epsilon)



    clf = CMC(clustering_clf, k=hyperParameters.k, m=hyperParameters.m)


    convoys = clf.fit_predict(data1)
    totalTime += time.time() - t1

    sorted_list = sorted(convoys, key=lambda x: x.start_time)
    printConvoys(sorted_list)


    t2 = time.time()
    atomTypes = calculateAtomTypes(fileLocation)
    finalConvoys = {}
    for index, convoy in enumerate(sorted_list):

        finalConvoys[index] = None

        arr = hydrogenBondCheckGivenIndices2(convoy.start_time, convoy.end_time, convoy.indices, data, atomTypes)

        finalConvoys[index] = arr


    totalTime += time.time() - t2



    #below this line we extract all BA convoys
    t3 = time.time()
    for i, indices in finalConvoys.items():
        print(i, consecutive_true_indices(indices))

    totalTime += time.time() - t3

    print(totalTime)


    #all final BAR convoys done below this line
    t4 = time.time()
    t = hyperParameters.t1
    data = consecutive_true_indices(finalConvoys[0])
    ranges = [(int(x.split('-')[0]), int(x.split('-')[1])) if '-' in x else (int(x), int(x)) for x in data.split(', ')]
    combined_ranges = combine_ranges(ranges, t + 1)

    # Output the result
    result = ""
    for i, (start, end) in enumerate(combined_ranges):
        result += f"{i+1}. {start}-{end}"
        if i < len(combined_ranges)-1:
            result += ", "
    print(result, time.time() - t4)

