import numpy as np
import numpy.linalg as la
from .Atom import Atom



class UnitCell:
    """Representation of a crystal unit cell."""
 
    atoms = []
    base = [[1,0,0],
            [0,1,0],
            [0,0,1]]
    rounding_decimals = 6

    def __init__(self, base=None, atoms=None):
        if base is None:
            base = np.array(self.__class__.base)
        self.base = base
        self.atoms = atoms or []        


    def _get_tofractionalcoordsMatrix(self):
        if not hasattr(self, '_tofractionalcoordsMatrix'):
            self._tofractionalcoordsMatrix = self._create_tofractionalcoordsMatrix()
        return self._tofractionalcoordsMatrix
    tofractionalcoords = property(_get_tofractionalcoordsMatrix)


    def _create_tofractionalcoordsMatrix(self):
        return la.inv(self.base.T)
    

    def calcFractionalCoords(self, v):
        '''compute the fractional coords of given vector'''
        tofractionalcoords = self.tofractionalcoords
        r = np.dot(tofractionalcoords, v)
        return np.round(r, self.rounding_decimals)


    def calcFractionalCoordsInCell(self, v):
        '''compute the fractional coords of given vector in the first cell'''
        v1 = self.calcFractionalCoords(v)
        return v1-np.floor(v1)


    def addAtom(self, atom):
        xyz = self.calcFractionalCoordsInCell(atom.xyz_cartn)
        symbol = atom.symbol
        atom = Atom(symbol, xyz)
        self.atoms.append(atom)
        return atom


    def hasAtom(self, atom):
        '''check if the given atom is alreday in the cell'''
        return self.getAtomIndex(atom) is not None


    def getAtomIndex(self, atom):
        '''get  the given atom is alreday in the cell'''
        xyz = atom.xyz_cartn
        xyz = self.calcFractionalCoordsInCell(xyz)
        for i, a in enumerate(self.atoms):
            if a.symbol == atom.symbol and (np.round(a.xyz - xyz, self.rounding_decimals)==[0,0,0]).all():
                return i
            continue
        return None


    def cartesian(self, u):
        """return cartesian coordinates of a lattice vector"""
        rc = np.dot(u, self.base)
        return rc


    def __str__(self):
        s_base = "base=%s" % self.base
        s_atoms = '\n'.join([str(a) for a in self.atoms])
        return s_base + '\n' + s_atoms



from numpy.testing import assert_array_almost_equal
import unittest
class TestCase(unittest.TestCase):

    def test1(self):
        base = np.array([[0,1,1], [1,0,1], [1,1,0]])
        uc = UnitCell(base=base)
        assert_array_almost_equal(
            uc.calcFractionalCoords([0,1,1]),
            [1,0,0])
        assert_array_almost_equal(
            uc.calcFractionalCoords([0,0.99999999,0.999999999]),
            [1,0,0])
        assert_array_almost_equal(
            uc.calcFractionalCoords([0,1.000000001,1.000000001]),
            [1,0,0])
        assert_array_almost_equal(
            uc.calcFractionalCoords([1,0,1]),
            [0,1,0])
        assert_array_almost_equal(
            uc.calcFractionalCoords([1,1,0]),
            [0,0,1])
        assert_array_almost_equal(
            uc.calcFractionalCoords([2,0,0]),
            [-1,1,1])
        
        assert_array_almost_equal(
            uc.calcFractionalCoordsInCell([0,1,1]),
            [0,0,0])
        assert_array_almost_equal(
            uc.calcFractionalCoordsInCell([1,0,1]),
            [0,0,0])
        assert_array_almost_equal(
            uc.calcFractionalCoordsInCell([1,1,0]),
            [0,0,0])
        assert_array_almost_equal(
            uc.calcFractionalCoordsInCell([2,0,0]),
            [0,0,0])
        return


    def test2(self):
        uc = UnitCell()
        
        uc.addAtom(Atom('H', (0,0,0.1)))
        
        self.assertTrue(uc.hasAtom(Atom('H', (0,0,0.10000001))))
        self.assertTrue(uc.hasAtom(Atom('H', (0,0,0.09999999))))
        self.assertTrue(not uc.hasAtom(Atom('H', (0,0,0.11))))
        self.assertTrue(not uc.hasAtom(Atom('H', (0,0,0.09))))

        self.assertTrue(uc.hasAtom(Atom('H', (0,0,1.10000001))))
        self.assertTrue(uc.hasAtom(Atom('H', (0,0,1.09999999))))
        self.assertTrue(not uc.hasAtom(Atom('H', (0,0,1.11))))
        self.assertTrue(not uc.hasAtom(Atom('H', (0,0,1.09))))
        self.assertTrue(not uc.hasAtom(Atom('H', (0,0,1.09999))))
        return



def main():
    unittest.main()
    return



if __name__ == '__main__': main()
