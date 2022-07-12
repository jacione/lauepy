# Olivier Delaire 

__doc__ = """Converters for UnitCell."""

from crystal.UnitCell import *
from crystal.Atom import *

def listOfAtom2UnitCell(loa):
    """Utility to convert a ListOfAtom instance to a UnitCell instance."""
    import ASE.ListOfAtoms
    import ASE.Atom
    
    symbols = [x.symbol for x in loa]
    cartpos = [x.position for x in loa]

    cellvectors = loa.cell
    uc = UnitCell()
    uc.setCellVectors(cellvectors)
    # note: by default, loa stores cartesian coordinates
    fracpos = [uc.cartesianToFractional(x) for x in cartpos]

    for n in range(len(symbols)):
        atom = Atom(symbol=symbols[n])
        uc.addAtom(atom, fracpos[n], '')

    return uc


def unitCell2ListOfAtom(uc, makeCartesian=False):
    """Utility to convert a UnitCell instance to a ASE ListOfAtom instance."""
    import ASE.ListOfAtoms
    import ASE.Atom
    loa = ASE.ListOfAtoms([],periodic=True) 
    loa.SetUnitCell(uc._cellvectors, fix=True)
    for site in uc:
        if makeCartesian:
            aseatom = ASE.Atom(site.getAtom().symbol, uc.fractionalToCartesian(site.getPosition().tolist()))
        else:
            aseatom = ASE.Atom(site.getAtom().symbol, site.getPosition().tolist())
        loa.append(aseatom)
    return loa

def p4vaspStruct2UnitCell(struct):
    """Utility to build a UnitCell instance from a  P4VASP structure.
    The atom types are not known from the POSCAR only (need POTCAR),
    but they are known if the Structure isntance is built from parsing
    the vasp.xml file properly, eg:
    from p4vasp.SystemPM import *
    from p4vasp.Structure import Structure
    run=XMLSystemPM('vasprun.xml')
    struct=run.FINAL_STRUCTURE
    """ 
    try:
        from AbInitio.vasp.parsing.Structure import Structure
    except ImportError:
        print("P4Vasp could not be imported in Python.")

    uc = UnitCell()
    cellvecs = [v.data for v in struct.basis]
    uc.setCellVectors(cellvecs)

    atomindex = 0
    for v in struct.positions:
        vaspatomstring = struct.getRecordForAtom(atomindex)['element']
        atomstring = vaspatomstring.split('_')[0]
        atom = Atom(symbol=atomstring.strip())
        site = Site(position=v.data, atom=atom)
        uc.addSite(site, siteId='')
        atomindex += 1
    return uc


def unitCell2P4vaspStruct(uc):
    """Create a P4VASP Structure instance from a UnitCell instance.
    The atom types and the number of atoms of each type need to be
    filled in from the UnitCell."""
    try:
        #from vaspparsing.SystemPM import *
        from AbInitio.vasp.parsing.Structure import Structure
        import AbInitio.vasp.parsing.matrix as p4mat
    except ImportError:
        print("P4Vasp could not be imported in Python.")
    struct = Structure()
    struct.basis = [p4mat.Vector(v.tolist()) for v in uc.getCellVectors()]
    #struct.positions = [p4mat.Vector(v.tolist()) for v in uc.getPositions()]
    # setup a dictionary of species index
    denum = uc.getAtomTypeDenum()
    atomtypes = list(denum.keys())
    indexdict = {}
    comment = ''
    for site in uc:
        atomtype = (site.getAtom().symbol, site.getAtom().mass)
        if atomtype not in indexdict:
            indexdict[atomtype] = atomtypes.index(atomtype)
            comment = comment + site.getAtom().symbol + str(denum[atomtype])
    struct.comment = comment
    for site in uc:
        pos = p4mat.Vector(site.getPosition().tolist())
        atomtype = (site.getAtom().symbol, site.getAtom().mass)
        speciesindex = indexdict[atomtype]
        #print atomtype, speciesindex
        struct.appendAtom(speciesindex, pos)
    return struct
    



    
