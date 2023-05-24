class ConvoyCandidate(object):
    """
    Attributes:
        indices(set): The object indices assigned to the convoy
        is_assigned (bool):
        start_time (int):  The start index of the convoy
        end_time (int):  The last index of the convoy
        atoms: these atoms represent the HB that is being tracked by this cluster/convoy. since we need to keep track of each HB and one cluster can have 2 different HB's. we need to keep them independently
    """
    __slots__ = ('indices', 'is_assigned', 'start_time', 'end_time', 'hb_hn', 'hb_n', 'hb_o', 'gap', 'totalGaps')

    def __init__(self, indices, is_assigned, start_time, end_time, hn = None, n = None, o = None, gap = 0, totalGaps = 0):
        self.indices = indices
        self.is_assigned = is_assigned
        self.start_time = start_time
        self.end_time = end_time
        self.hb_hn = hn 
        self.hb_n = n
        self.hb_o = o
        self.gap = gap
        self.totalGaps = totalGaps

    def __repr__(self):
        return '<%r %r indices=%r, is_assigned=%r, start_time=%r, end_time=%r, hb_hn=%r, hb_n=%r, hb_o=%r, gap=%r, totalGaps=%r>' % (self.__class__.__name__, id(self), self.indices, self.is_assigned, self.start_time, self.end_time, self.hb_hn, self.hb_n, self.hb_o, self.gap, self.totalGaps)
