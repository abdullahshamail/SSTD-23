class Cluster(object):
    """
Attributes:
    cluster_id: the cluster number this cluster was assigned by the algorithm
    timestep: the timestep this cluster belongs to
    hb_present: does this cluster have a hydroegn bond present or not
    hb_list: the list of hydrogen bonds
"""
    __slots__ = ('cluster_id', 'timestep', 'hb_present', 'hb_list', 'indices', 'totalHB')

    def __init__(self, cluster_id, timestep, hb_present, hb_list, indices, totalHB = 0):
        self.cluster_id = cluster_id
        self.timestep = timestep
        self.hb_present = hb_present
        self.hb_list = hb_list
        self.indices = indices
        self.totalHB = totalHB

    def __repr__(self):
        return '<%r %r cluster_id=%r, timestep=%r, hb_present=%r, hb_list=%r, indices=%r, totalHB=%r>' % (self.__class__.__name__, id(self), self.cluster_id, self.timestep, self.hb_present, self.hb_list, self.indices, self.totalHB)