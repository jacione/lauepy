
'''
find out the Strukturbericht Designation of a structure

now only handles bcc and fcc

TODO: add much more handlers. 
'''


class StrukturberichtDesignation:

    def __init__(self, symbol = None, alias = None):
        self.symbol = symbol
        self.alias = alias
        return


    def __str__(self):
        return '%s(%s)' % (self.symbol, self.alias)
    __repr__ = __str__


    def __eq__(self, rhs):
        return self.symbol == rhs.symbol



A1 = StrukturberichtDesignation(symbol='A1', alias='fcc')
A2 = StrukturberichtDesignation(symbol='A2', alias='bcc')
Ah = StrukturberichtDesignation(symbol='Ah', alias='simple cubic')
B1 = StrukturberichtDesignation(symbol='B1', alias='NaCl')
B2 = StrukturberichtDesignation(symbol='B2', alias='CsCl')


class StrukturberichtDesignationFinder:

    def find(self, struct):
        n = len(struct)

        sg = struct.sg
        number = sg.number
        
        handler = '_on_%s_%s' % (number, n)
        
        if not hasattr(self, handler): return
        method = getattr(self, handler)
        return method(struct)


    def _on_221_1(self, struct):
        return Ah


    def _on_221_2(self, struct):
        lattice = struct.lattice
        if lattice.alpha == 90 and lattice.beta == 90 and lattice.gamma == 90 \
               and struct[0].symbol != struct[1].symbol:
            return B2


    def _on_229_1(self, struct):
        lattice = struct.lattice
        if isAlmostEqual(lattice.alpha, acos1_3) and \
           isAlmostEqual(lattice.beta, acos1_3) and \
           isAlmostEqual(lattice.gamma, acos1_3):
            return A2

    
    def _on_229_2(self, struct):
        lattice = struct.lattice
        if lattice.alpha == 90 and lattice.beta == 90 and lattice.gamma == 90 \
               and struct[0].symbol == struct[1].symbol:
            return A2
        return


    def _on_225_1(self, struct):
        lattice = struct.lattice
        if lattice.alpha == 60 and lattice.beta == 60 and lattice.gamma == 60:
            return A1


    def _on_225_2(self, struct):
        lattice = struct.lattice
        if lattice.alpha == 60 and lattice.beta == 60 and lattice.gamma == 60 \
            and struct[0].symbol != struct[1].symbol:
            return B1


def isAlmostEqual(a, b, places=5):
    return round(abs(b-a), places) == 0


import math
acos1_3 = math.acos(1./3)*180/math.pi



from danse.ins import matter

import unittest
class TestCase(unittest.TestCase):

    def testA1(self):
        struct = matter.Structure(
            lattice = matter.Lattice(a=1,b=1,c=1,alpha=60,beta=60,gamma=60),
            sgid = 225,
            atoms = [matter.Atom('Cu')],
            )
        self.assertEqual(StrukturberichtDesignationFinder().find(struct), A1)
        return


    def testA2(self):
        struct = matter.Structure(
            lattice = matter.Lattice(a=1,b=1,c=1,alpha=acos1_3,beta=acos1_3,gamma=acos1_3),
            sgid = 229,
            atoms = [matter.Atom('Cu')],
            )
        self.assertEqual(StrukturberichtDesignationFinder().find(struct), A2)

        struct = matter.Structure(
            lattice = matter.Lattice(a=1,b=1,c=1,alpha=90,beta=90,gamma=90),
            sgid = 229,
            atoms = [matter.Atom('Cu', xyz=[0,0,0]), matter.Atom('Cu', xyz=[0.5,0.5,0.5])],
            )
        self.assertEqual(StrukturberichtDesignationFinder().find(struct), A2)

        return


    def testAh(self):
        struct = matter.Structure(
            lattice = matter.Lattice(a=1,b=1,c=1,alpha=90,beta=90,gamma=90),
            sgid = 221,
            atoms = [matter.Atom('Cu')],
            )
        self.assertEqual(StrukturberichtDesignationFinder().find(struct), Ah)
        return


    def testB1(self):
        struct = matter.Structure(
            lattice = matter.Lattice(a=1,b=1,c=1,alpha=60,beta=60,gamma=60),
            sgid = 225,
            atoms = [matter.Atom('Na'), matter.Atom('Cl', xyz=[0.5,0.5,0.5])],
            )
        self.assertEqual(StrukturberichtDesignationFinder().find(struct), B1)
        return


    def testB2(self):
        struct = matter.Structure(
            lattice = matter.Lattice(a=1,b=1,c=1,alpha=90,beta=90,gamma=90),
            sgid = 221,
            atoms = [matter.Atom('Cs'), matter.Atom('Cl', xyz=[0.5,0.5,0.5])],
            )
        self.assertEqual(StrukturberichtDesignationFinder().find(struct), B2)
        return


    def testMethodInStructureClass(self):
        struct = matter.Structure(
            lattice = matter.Lattice(a=1,b=1,c=1,alpha=60,beta=60,gamma=60),
            sgid = 225,
            atoms = [matter.Atom('Cu')],
            )
        self.assertEqual(struct.StrukturberichtDesignation, A1)
        return



def main():
    unittest.main()
    return



if __name__ == '__main__': main()
