# this class does not really seem to be necessary with lattice
class Site:
    """Representation of a crystallographic site.
    A site corresponds to a position in fractional coordinates.
    It also may have an atom associated with it."""

    def __init__(self, position='', atom=None, occproba=1.0):

        for x in position:
            if (x < 0) or (x > 1):
                raise ValueError("Site coordinates must be fractional positions.")
        self._position = np.array( position )
        self._atom = atom
        self._occproba = occproba
        self.xyz=self._position

    def __str__(self):
        rt = str(self._position) + ":" + str(self._atom)
        return rt

    def setPosition(self, position):
        """Sets the position (in fractional coordinates) of a site.""" 
        for x in position:
            if (x < 0) or (x > 1):
                raise ValueError("Site coordinates must be fractional positions.")
        self._position = np.array(position)
        
    def getPosition(self):
        """Get the position (in fractional coordinates) of a site."""
        return self._position

    def setAtom(self, atom):
        """Set an atom at a site."""
        self._atom = atom

    def getAtom(self):
        """Returns the atom at a site."""
        return self._atom

    def getOccProba(self):
        """Returns the occupation probability for the atom at this site."""
        return self._occproba

    def setOccProba(self,p):
        """Sets the occupation probability for the atom at this site."""
        try:
            proba = float(p)
        except:
            raise ValueError('Probability should be a number in [0,1].')
        if proba <= 1.0 and proba >= 0.0:
            self.__occproba = proba
        else:
            raise ValueError('Probability should be a number in [0,1].')
            
        
    pass # end of class site




