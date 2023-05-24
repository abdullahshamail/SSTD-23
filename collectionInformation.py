from cluster import Cluster
from utilities import getHBinformation
import numpy as np
from convoycandidate import ConvoyCandidate
class CollectionInformation(object):
    """
    Attributes:
        column(set): the column number (timestep) that this object contains clusters for
        clusters: clusters for the timestep
        unique-clusters:  self explanatory
        cluster_indices:  ConvoyCandidate object that contains the cluster atom indices
    """
    __slots__ = ('column', 'clusters', 'unique_clusters', 'clusters_indices', 'hydrogen_bonds', 'num_bonds', 'HBClusters')

    def __init__(self, column, clusters, unique_clusters, clusters_indices):
        self.column = column
        self.clusters = clusters
        self.unique_clusters = unique_clusters
        self.clusters_indices = clusters_indices #a dictinoary of ConvoyCandidate class
        self.hydrogen_bonds = [] #each has a convoyCandidate object
        self.HBClusters = []
        self.num_bonds = 0 #default is zero, total number of bonds in the collection

    def __repr__(self):
        return '<%r %r column=%r, clusters=%r, unique_clusters=%r, clusters_indices=%r, hydrogen_bonds=%r>' % (self.__class__.__name__, id(self), self.column, self.clusters, self.unique_clusters, self.clusters_indices, self.hydrogen_bonds)
    
    def calculateHB(self, data):
        for i, candidate in self.clusters_indices.items():
            temporary_cluster = Cluster(cluster_id = i, timestep = self.column, hb_present=False, hb_list=[], indices=np.array(list(candidate.indices)), totalHB=0)
            updatedCluster = getHBinformation(data, temporary_cluster, self.column)
            convoy_cand = ConvoyCandidate(indices=candidate.indices, is_assigned=False, start_time=None, end_time=None)
            if updatedCluster.hb_present:
                self.hydrogen_bonds.append(convoy_cand)
                self.HBClusters.append(updatedCluster)
            self.num_bonds += updatedCluster.totalHB