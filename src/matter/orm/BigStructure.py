# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                   Jiao Lin
#                      California Institute of Technology
#                      (C) 2006-2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from dsaw.model.Inventory import Inventory as InvBase


# data object
from ..Structure import Structure


# inventory
from dsaw.model.Inventory import Inventory as InvBase
class Inventory(InvBase):

    id = InvBase.d.str(name="id", max_length=64, constraints = 'PRIMARY KEY')
    short_description = InvBase.d.str(
        name = 'short_description', max_length = 256, default ="", label='Description')
    lattice_base = InvBase.d.array(name='lattice_base', shape=(3,3), elementtype='float')
    atom_symbols = InvBase.d.array(name='atom_symbols', elementtype='str')
    atom_positions = InvBase.d.array(name='atom_positions', elementtype='float')
    spacegroup_num = InvBase.d.int(name = 'spacegroup_num', default =1, label='Spacegroup #')
    chemical_formula = InvBase.d.str(name='chemical_formula', max_length=1024)
    primitive_cell_base = InvBase.d.array(name='primitive_cell_base', shape=(3,3), elementtype='float')
    date = InvBase.d.date(name='date')

    dbtablename = 'bigstructure'

Structure.Inventory = Inventory


# dsaw.model helpers
def __establishInventory__(self, inventory):
    inventory.short_description = self.description
    inventory.lattice_base = self.lattice.base
    inventory.atom_symbols = self.symbols
    inventory.atom_positions = self.xyz
    inventory.spacegroup_num = self.sg.number
    inventory.chemical_formula = self.getChemicalFormula()
    inventory.primitive_cell_base = self.primitive_unitcell.base
    return
Structure.__establishInventory__ = __establishInventory__

def __restoreFromInventory__(self, inventory):
    atoms=[]
    from matter.Atom import Atom
    from matter.Lattice import Lattice
    for symbol,pos in zip(inventory.atom_symbols, inventory.atom_positions):
        atoms.append(Atom(symbol, pos))
    self.__init__(atoms = atoms,
                  lattice = Lattice(base = inventory.lattice_base),
                  sgid = inventory.spacegroup_num,
                  description = inventory.short_description,
                  )
    #self.primitive_unitcell = inventory.primitive_unitcell
    return
Structure.__restoreFromInventory__ = __restoreFromInventory__


# version
__id__ = "$Id$"

# End of file 
