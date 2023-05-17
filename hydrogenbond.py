class HydrogenBond(object):
    """
    Attributes:
        hydrogen_index: index of hydrogen
        nitrogen_index: index of nitrogen
        oxygen_index: index of oxygen
        timeframe: timeframe for hbond
        *Coord = coordinates of the respective atom at that timestep
        oxygenType = type of oxygen
    """
    __slots__ = ('hydrogen_index', 'nitrogen_index', 'oxygen_index', 'timeframe', 'hydrogenCoord', 'nitrogenCoord', 'oxygenCoord', 'oxygenType')

    def __init__(self, hydrogen_index, hydrogenCoord, nitrogen_index, nitrogenCoord, oxygen_index, oxygenCoord, timeframe, oxgyenType):
        self.hydrogen_index = hydrogen_index
        self.nitrogen_index = nitrogen_index
        self.oxygen_index = oxygen_index
        self.timeframe = timeframe
        self.hydrogenCoord = hydrogenCoord
        self.nitrogenCoord = nitrogenCoord
        self.oxygenCoord = oxygenCoord
        self.oxygenType = oxgyenType

    def __repr__(self):
        return '<%r %r hydrogen_index=%r, hydrogenCoord=%r, nitrogen_index=%r, nitrogenCoord=%r, oxygen_index=%r, oxygenCoord=%r, timeframe=%r, oxygenType=%r>' % (self.__class__.__name__, id(self), self.hydrogen_index, self.hydrogenCoord, self.nitrogen_index, self.nitrogenCoord, self.oxygen_index, self.oxygenCoord, self.timeframe, self.oxygenType)
    
    def hb_atom_indexes(self):
        return self.hydrogen_index, self.nitrogen_index, self.oxygen_index