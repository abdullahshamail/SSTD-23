# Import all libraries

import numpy as np
from sklearn.cluster import DBSCAN
import pandas as pd
from convoycandidate import ConvoyCandidate
from collectionInformation import CollectionInformation
from utilities import HB_truth, getCoordinates, find_key_by_value, hydrogenBondCheckGivenIndices, calculateAtomTypes, filterDF, extractInterestingMoleculesPerTimeframe, printConvoys
from hyperparameters import Hyperparameters
import time



class SemanticsBAR(object):
    """semantic aware algorithm for detecting bar-convoys

    Attributes:
        k (int):  Min number of consecutive timestamps to be considered a convoy
        m (int):  Min number of elements to be considered a convoy
        t1 (int): Max allowed time gap
        t2 (int): total allowed time gaps
    """
    def __init__(self, clf, k, m, t1, t2, data, threshold, df, atomType, nitrogens_and_oxygens_indices):
        self.clf = clf
        self.k = k
        self.m = m
        self.t1 = t1
        self.t2 = t2
        self.data = data
        self.threshold = threshold
        self.df = df
        self.atomType = atomType
        self.nitrogens_and_oxygens_indices = nitrogens_and_oxygens_indices

    def HBCheck(self, atomsIHave, column):
        for hn in atomsIHave["hn"]:
            for n in atomsIHave["n"]:
                coord_hn = getCoordinates(self.data[column], hn)
                coord_n = getCoordinates(self.data[column], n)
                for o in atomsIHave["o"]:
                    coord_o = getCoordinates(self.data[column], o)
                    if HB_truth(coord_hn, coord_n, coord_o):
                        return True
        
        return False

    def createClusters(self, X, y, column, sample_weight, indicesOfFilteredData):
        values = [row[column] if isinstance(row[column], (list, set)) else [row[column]] for row in X]
        values = np.array(values)[indicesOfFilteredData]
        if len(values) < self.m:
            return 1, 0 , 0 , 0
        clusters = self.clf.fit_predict(values, y=y, sample_weight=sample_weight)
        unique_clusters = set(clusters)
        clusters_indices = dict((cluster, ConvoyCandidate(indices=set(), is_assigned=False, start_time=None, end_time=None)) for cluster in unique_clusters)

        return 0, clusters, unique_clusters, clusters_indices
         
    def fit_predict(self, X, y=None, sample_weight=None):
        convoy_candidates = set()
        columns = len(X[0])
        column_iterator = range(columns)
        output_convoys = []

        for column in column_iterator:

            results = extractInterestingMoleculesPerTimeframe(self.data[column], self.nitrogens_and_oxygens_indices, self.threshold)
            reducedMolecules = []
            if len(results) > 0:
                reducedMolecules = np.unique(np.hstack((np.unique(np.array(results)[:, 0]), np.unique(np.array(results)[:, 1]))))

            if len(reducedMolecules) > 0:
                indicesOfFilteredData = np.sort(self.df[self.df["subst_id"].isin(reducedMolecules)].reset_index()["index"].values)
                tempIndices = np.arange(0, indicesOfFilteredData.shape[0])
                reverser = dict(zip(tempIndices, indicesOfFilteredData)) #reverser to use to get original indices of the atoms in clusters after cluster creation

            current_convoy_candidates = set()

            enough_objects, clusters, unique_clusters, clusters_indices = self.createClusters(X, y, column, sample_weight, indicesOfFilteredData)

            if enough_objects == 1 or len(reducedMolecules) == 0:
                continue


            forwardClusters = {}
            for index, cluster_assignment in enumerate(clusters):
                clusters_indices[cluster_assignment].indices.add(reverser[index])

            forwardClusters[column] = CollectionInformation(column = column, clusters = clusters, unique_clusters = unique_clusters, clusters_indices = clusters_indices)

            for forward_col in range(column + 1, column + 1 + self.t1):
                if forward_col < columns:
                    results_f = extractInterestingMoleculesPerTimeframe(self.data[forward_col], self.nitrogens_and_oxygens_indices, self.threshold)
                    reducedMolecules_f = []
                    if len(results_f) > 0:
                        reducedMolecules_f = np.unique(np.hstack((np.unique(np.array(results_f)[:, 0]), np.unique(np.array(results_f)[:, 1]))))

                    if len(reducedMolecules_f) > 0:
                        indicesOfFilteredData_f = np.sort(self.df[self.df["subst_id"].isin(reducedMolecules_f)].reset_index()["index"].values)
                        tempIndices_f = np.arange(0, indicesOfFilteredData_f.shape[0])
                        reverser_f = dict(zip(tempIndices_f, indicesOfFilteredData_f))
                        
                    enough_objects_f, clusters_f, unique_clusters_f, clusters_indices_f = self.createClusters(X, y, forward_col, sample_weight, indicesOfFilteredData_f)


                    for index_f, cluster_assignment_f in enumerate(clusters_f):

                        clusters_indices_f[cluster_assignment_f].indices.add(reverser_f[index_f])

                    forwardClusters[forward_col] = CollectionInformation(column = forward_col, clusters = clusters_f, unique_clusters = unique_clusters_f, clusters_indices = clusters_indices_f)


            # update existing convoys
            for convoy_candidate in convoy_candidates:

                found = False
                for curr_col in range(column, column + 1 + self.t1):
                    if curr_col < columns and not found:

                        convoy_candidate.is_assigned = False
                        for cluster in forwardClusters[curr_col].unique_clusters:
                            cluster_indices = forwardClusters[curr_col].clusters_indices[cluster].indices
                            cluster_candidate_intersection = cluster_indices & convoy_candidate.indices
                            if len(cluster_candidate_intersection) < self.m:
                                continue

                            atomsIHave = {
                                "n": [],
                                "hn": [],
                                "o": []
                            }
                            for index in cluster_candidate_intersection:
                                if not find_key_by_value(index, self.atomType) == None:
                                    atomsIHave[find_key_by_value(index, self.atomType)].append(index)
                            if self.HBCheck(atomsIHave, curr_col):
                                convoy_candidate.indices = cluster_candidate_intersection
                                # print("redo",convoy_candidate.indices)
                                convoy_candidate.end_time = curr_col
                                forwardClusters[curr_col].clusters_indices[cluster].is_assigned = convoy_candidate.is_assigned = True
                                found = True
                                current_convoy_candidates.add(convoy_candidate)
                        if not found:
                            convoy_candidate.totalGaps += 1

                # check if candidates qualify as convoys
                candidate_life_time = (convoy_candidate.end_time - convoy_candidate.start_time) + 1

                if (not convoy_candidate.is_assigned or column == column_iterator[-1]) and candidate_life_time >= self.k and convoy_candidate.totalGaps <= self.t2:

                    output_convoys.append(convoy_candidate)

            # create new candidates
            for cluster in unique_clusters:
                cluster_data = clusters_indices[cluster]

                atomsIHave = {
                    "n": [],
                    "hn": [],
                    "o": []
                }
                for index in cluster_data.indices:
                    if not find_key_by_value(index, self.atomType) == None:
                        atomsIHave[find_key_by_value(index, self.atomType)].append(index)
                if cluster_data.is_assigned:
                    continue
                if not self.HBCheck(atomsIHave, column):
                    continue
                
                cluster_data.start_time = cluster_data.end_time = column
                current_convoy_candidates.add(cluster_data)
            convoy_candidates = current_convoy_candidates
        return output_convoys
    

hyperParameters = Hyperparameters(m = 50, k = 5, epsilon = 3.5, t1 = 1, t2 = 50, threshold=4)
dataSize = 500

if __name__ == '__main__':
    

    data = np.load(r"..\xyzCoords.npy")

    df = pd.read_csv(r"...\molecularInformation.csv") #location for information of atoms

    nitrogens_and_oxygens_indices = filterDF(df)


    # Transpose the data for Convoy algorithm

    hyperParameters = Hyperparameters(m = 50, k = 10, epsilon = 3.5, t = 1, threshold=4)
    transposedData = list()
    for x in range(data.shape[1]):
        transposedData.append(data[0:dataSize,x,0:3].tolist())

    np.array(transposedData).shape


    # get indices of atoms with atom_type 'n' and subst_id between 1-4
    n_atoms = df[(df['atom_type'] == 'n') & (df['subst_id'] >= 1) & (df['subst_id'] <= 4)].index.tolist()
    hn_atoms = df[(df['atom_type'] == 'hn') & (df['subst_id'] >= 1) & (df['subst_id'] <= 4)].index.tolist()

    # get indices of atoms with atom_type 'o' or 'os' and subst_id between 5-14
    o_atoms = df[((df['atom_type'] == 'o') | (df['atom_type'] == 'os')) & (df['subst_id'] >= 5) & (df['subst_id'] <= 14)].index.tolist()

    atomType = {
        "n": n_atoms,
        "hn": hn_atoms,
        "o": o_atoms
    }

    totalTime = 0
    t1 = time.time()
    clustering_clf = DBSCAN(eps=hyperParameters.epsilon)

    clf = SemanticsBAR(clustering_clf, k=hyperParameters.k, m=hyperParameters.m, t1=hyperParameters.t, t2=hyperParameters.t2, data=data, threshold=hyperParameters.threshold, df=df, atomType=atomType, nitrogens_and_oxygens_indices=nitrogens_and_oxygens_indices)


    convoys = clf.fit_predict(transposedData)

    totalTime += time.time() - t1
    sorted_list = sorted(convoys, key=lambda x: x.start_time)

    print(totalTime)
    printConvoys(sorted_list)
