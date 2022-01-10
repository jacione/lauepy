# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""ChemicalElements module

This module defines a function `Element()`.  The argument to the
`Element()` function can be a chemical symbol or an atomic number.
The object returned has the following attributes: ``symbol``, ``number``,
``name``, and if it makes sence also ``mass``, ``covalent_radius`` and
``crystal_structure``.

>>> m = Element('Si').mass
>>> r = Element('Fe').covalent_radius

"""

__docformat__ = 'reStructuredText'



#from ASE import enumerate
#from ASE.ChemicalElements.symbol import symbols



symbols = ['X',  'H',  'He', 'Li', 'Be',
           'B',  'C',  'N',  'O',  'F',
           'Ne', 'Na', 'Mg', 'Al', 'Si',
           'P',  'S',  'Cl', 'Ar', 'K',
           'Ca', 'Sc', 'Ti', 'V',  'Cr',
           'Mn', 'Fe', 'Co', 'Ni', 'Cu',
           'Zn', 'Ga', 'Ge', 'As', 'Se',
           'Br', 'Kr', 'Rb', 'Sr', 'Y',
           'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
           'Rh', 'Pd', 'Ag', 'Cd', 'In',
           'Sn', 'Sb', 'Te', 'I',  'Xe',
           'Cs', 'Ba', 'La', 'Ce', 'Pr',
           'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
           'Tb', 'Dy', 'Ho', 'Er', 'Tm',
           'Yb', 'Lu', 'Hf', 'Ta', 'W',
           'Re', 'Os', 'Ir', 'Pt', 'Au',
           'Hg', 'Tl', 'Pb', 'Bi', 'Po',
           'At', 'Rn', 'Fr', 'Ra', 'Ac',
           'Th', 'Pa', 'U',  'Np', 'Pu',
           'Am', 'Cm', 'Bk', 'Cf', 'Es',
           'Fm', 'Md', 'No', 'Lw']


numbers = {}
"""A dictionary translating chemical symbols to atomic numbers."""

# Fill in the numbers:
for _Z, _symbol in enumerate(symbols):
    numbers[_symbol] = _Z

_elements = {}


def Element(Z):
    if type(Z) is str:
        Z = numbers[Z]
    try:
        return _elements[Z]
    except KeyError:
        element = _Element(Z)
        _elements[Z] = element
        return element
    

class _Element:
    def __init__(self, Z):
        self.number = Z
        
    def __getattr__(self, name):
        try:
            module = __import__(name, globals(), locals(), [])
        except ImportError:
            raise AttributeError('Unknown property!')
        stuff = getattr(module, '_data')[self.number]
        if stuff is None:
            raise RuntimeError('Unknown!')
        setattr(self, name, stuff)
        return stuff
