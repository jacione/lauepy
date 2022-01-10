########################################################################
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

import numpy
from .Atom import Atom
from .Lattice import Lattice

class Structure(list):
    """Structure --> group of atoms

    Structure class is inherited from Python list.  It contains
    a list of Atom instances.  Structure overloads setitem and setslice
    methods so that the lattice attribute of atoms get set to lattice.

    Data members:
        description   -- structure description
        lattice -- coordinate system (instance of Lattice)
    """
    

    def __init__(self, atoms=[], lattice=None, sgid=1, description=None, filename=None,
                 primitive_unitcell=None
                 ):
        """define group of atoms in a specified lattice.

        atoms    -- list of Atom instances to be included in this Structure.
                    When atoms argument is an existing Structure instance,
                    the new Structure is its copy.
        lattice  -- instance of Lattice defining coordinate systems, property.
        sgid     -- space group symbol, either short_name or pdb_name,
                    whatever it means in mmlib.  Can be also an integer.
        description -- string description of the structure
        filename -- optional, name of a file to load the structure from.
                    Overrides atoms argument when specified.

        Structure(stru)     create a copy of Structure instance stru.

        Because Structure is inherited from a list it can use list expansions,
        for example:
            oxygen_atoms = [ for a in stru if a.symbol == "O" ]
            oxygen_stru = Structure(oxygen_atoms, lattice=stru.lattice)
        """

        self._labels = {}
        self._labels_cached = False
        if isinstance(atoms, Structure):
            stru = atoms
            # create a shallow copy of all source attributes
            self.__dict__.update(stru.__dict__)
            # make a deep copy of source lattice
            self.lattice = Lattice(stru.lattice)

        self.description = description
        self.sg = sgid
        # check if data should be loaded from file
        if filename is not None:
            self.read(filename)
        # otherwise assign list of atoms to self
        else:
            self[:] = atoms
            
        # override from lattice argument
        if lattice is None:
            if not self.lattice:    self.lattice = Lattice()
        #elif not isinstance(lattice, Lattice):
        #   emsg = "expected instance of Lattice"
        #    raise TypeError(emsg)
        else:
            self.lattice = lattice


        self._primitive_unitcell = primitive_unitcell
        import time
        self.date = time.ctime()


    def __str__(self):
        """simple string representation"""
        s_lattice = "lattice=%s" % self.lattice
        s_atoms = '\n'.join([str(a) for a in self])
        return s_lattice + '\n' + s_atoms


    def getChemicalFormula(self):
        atoms = self
        counts = {}
        for atom in atoms:
            e = atom.symbol
            if e in counts: counts[e]+=1
            else: counts[e]=1
            continue
        elems = list(counts.keys())
        elems.sort()
        chemFormRaw = ''.join( '%s_%s ' % (e, counts[e]) for e in elems )
        return chemFormRaw.strip()
    
    def getSpecies(self):
        speciesList = []
        for atom in self:
            if atom.symbol in speciesList:
                pass
            else:
                speciesList.append(atom.symbol)
        return speciesList

    # these properties are actually independent of space group and so should
    # be moved to lattice perhaps
    def _getCrystalSystem(self):
        return self.sg.crystal_system
    crystal_system = property(_getCrystalSystem)


    def _getCentering(self):
        return self.sg.pdb_name[0]
    centering = property(_getCentering)

    
    _centerings = {
        'P': 'primitive',
        'C': 'single face centered',
        'A': 'single face centered',
        'H': 'primitive',
        'R': 'primitive',
        'I': 'body centered',
        'F': 'face centered',
        }
    def _getCenteringDescription(self):
        return self._centerings[self.centering]
    centering_description = property(_getCenteringDescription)
    
    def _getBravaisType(self):
        centeringAndSystem2bravais = {
        ('P', 'CUBIC'):'simple cubic',
        ('F', 'CUBIC'):'face centered cubic',
        ('I', 'CUBIC'):'body centered cubic', 
        ('P', 'HEXAGONAL'):'simple hexagonal',
        ('P', 'TRIGONAL'):'simple trigonal',
        ('P','TETRAGONAL'):'simple tetragonal',
        ('R','TETRAGONAL'):'simple tetragonal',
        ('I','TETRAGONAL'):'body centered tetragonal',
        ('P','ORTHORHOMBIC'):'simple orthorhombic',
        ('C','ORTHORHOMBIC'):'base centered orthorhombic',
        ('A','ORTHORHOMBIC'):'base centered orthorhombic',
        ('F','ORTHORHOMBIC'):'face centered orthorhombic',
        ('I','ORTHORHOMBIC'):'body centered orthorhombic',
        ('P','ORTHORHOMBIC'):'simple monoclinic',
        ('P','MONOCLINIC'):'base centered monoclinic',
        ('P','TRICLINIC'):'triclinic'                  
        }
        return centeringAndSystem2bravais[(self.sg.short_name[0], self.sg.crystal_system)]
    bravais_type = property(_getBravaisType)


    def _getStrukturberichtDesignation(self):
        finder = self._getStrukturberichtDesignationFinder()
        return finder.find(self)
    StrukturberichtDesignation = property(_getStrukturberichtDesignation)
        

    def _getStrukturberichtDesignationFinder(self):
        k = '_StrukturberichtDesignationFinder'
        if hasattr(self, k): return getattr(self, k)
        from .StrukturberichtDesignationFinder import StrukturberichtDesignationFinder
        r = StrukturberichtDesignationFinder()
        setattr(self, k, r)
        return r
    

    def _get_primitive_unitcell(self):
        t = self._primitive_unitcell
        if t is None:
            #t = self._primitive_unitcell = self._lattice.getPrimitiveLattice(self.sg)
            t = self._primitive_unitcell = self._create_primitive_unitcell()
        return t
    def _set_primitive_unitcell(self, p):
        self._primitive_unitcell = p
        return p
    primitive_unitcell = property(_get_primitive_unitcell, _set_primitive_unitcell)


    def _create_primitive_unitcell(self):
        # the ctor
        from .UnitCell import UnitCell
        base = self._lattice.getPrimitiveLattice(self.sg)
        base = numpy.array(base)
        
        puc = UnitCell(base=base)
        for atom in self:
            if puc.hasAtom(atom): continue
            puc.addAtom(atom)
            continue                
        return puc
        
        # check symmetry
        sg = self.sg
        verdict, atompos, symop = self.symConsistent()
        if not verdict:
            raise RuntimeError("inconsistent structure: atom position: %s, sym operation %s" % (
                atompos, symop))

        #
        if sg.number == 225:
            # face centered--below is only for cubic...need to generalize
            a = self.lattice.a
            base = numpy.array([[0,1,1], [1,0,1], [1,1,0]])*0.5*a
            
        elif sg.number == 229:
            # body centered--below is only for cubic...need to generalize
            a = self.lattice.a
            base = numpy.array([[1,0,0], [0,1,0], [0.5,0.5,0.5]])*a
            
            
        elif sg.short_name[0] == 'P':
            # primitive
            base = self.lattice

        else:
            # not implemented
            return
        
        
    
    def addNewAtom(self, *args, **kwargs):
        """Add new Atom instance to the end of this Structure.

        All arguments are forwarded to Atom constructor.
        """
        kwargs['lattice'] = self.lattice
        a = Atom(*args, **kwargs)
        list.append(self, a)
        self._uncache('labels')

    def getLastAtom(self):
        """Return Reference to the last Atom in this structure.
        """
        last_atom = self[-1]
        return last_atom

    def getAtom(self, id):
        """Reference to internal Atom specified by the identifier.

        id  -- zero based index or a string label formatted as
               "%(element)s%(order)i", for example "Na1", "Cl1"

        Return Atom instance.
        Raise ValueError for invalid id.

        See also getLabels().
        """
        try:
            if type(id) is int:
                rv = self[id]
            else:
                if not self._labels_cached or id not in self._labels:
                    self._update_labels()
                rv = self._labels[id]
        except (IndexError, KeyError):
            emsg = "Invalid atom identifier %r." % id
            raise ValueError(emsg)
        return rv


    def getLabels(self):
        """List of unique string labels for all atoms in this structure.

        Return a list.
        """
        elnum = {}
        labels = []
        for a in self:
            elnum[a.symbol] = elnum.get(a.symbol, 0) + 1
            alabel = a.symbol + str(elnum[a.symbol])
            labels.append(alabel)
        return labels
        
    #untested--from UnitCell
    def bringFractionalPositionIntoCell(self, fracpos):
        """Brings a fractional position (x,y,z) 'into' the unit cell,
        i.e.: (x,y,z)->(x',y',z') such that x,y,z in [0,1( """
        pos = numpy.array(fracpos)
        assert (len(pos) == 3)
        for i in range(3):
            if pos[i]<0:
                while pos[i]<0:
                    pos[i] += 1
            if pos[i]>=1:
                while pos[i]>=1:
                    pos[i] -= 1
        return pos

    #untested--from UnitCell
    def cartesianPositionInLattice(self, siteId, latticeVector):
        """Returns the Cartesian position vector from the origin
        ( fractional coordinates [0,0,0] in unit cell [0,0,0]),
        for a Site corresponding to 'siteId',
        in the unit cell corresponding to latticeVector
        (triplets of coordinates in terms of cellvectors), 
        defining which unit cell in the lattice.
        """
        try:
            posincell = self.getCartesianPosition(siteId)
        except KeyError: raise KeyError('Invalid site Id')
        pos = numpy.array(posincell) + numpy.dot(latticeVector, self._lattice)
        return pos
    
    def getPosition(self, siteId):
        """Returns the (fractional) position of a site."""
        return self._siteIds[siteId].getPosition()

    def setPositions(self, positions):
        """Sets the (fractional) positions of the sites in the unit cell."""
        assert(len(positions) == self.getNumSites())
        for isite in range(self.getNumSites()):
            self._sites[isite].setPosition(positions[isite])

    def getCartesianPosition(self, siteId):
        """Returns the cartesian position of a site."""
        return self.fractionalToCartesian(self._siteIds[siteId].getPosition())

    def fractionalToCartesian(self, fracpos):
        """Converts a coordinate from fractional to cartesian.""" 
        return (fracpos * self._lattice).sum(0)  # should be double-checked

    def cartesianToFractional(self, cartpos):
        """Converts a coordinate from cartesian to fractional."""
        return (cartpos * numpy.linalg.inv(self._lattice)).sum(0)  # should be double-checked
    
    def generateDescription(self):
        if self._description==None:
            self._description = self.getChemicalFormula()#+' in '+str(self.lattice)
        return self._description
    def setDescription(self, desc):
        self._description = desc
    description = property(generateDescription, setDescription, "structure description")
    
################################################    
# property methods
################################################
#Notes:
# for now these are done in the style of diffraction
# eventually will be done in Jiao's style with metaclasses

    # fractional xyz
    def _get_xyz(self):
        return [atom.xyz.tolist()  for atom in self[:]]
    def _set_xyz(self, xyzList):
        for atom,xyz in zip(self, xyzList):
            atom.xyz = xyz
    xyz = property(_get_xyz, _set_xyz, doc =
        """fractional coordinates of all atoms""" )  

    # xyz_cartn
    def _get_xyz_cartn(self):
        return [atom.xyz_cartn.tolist() for atom in self[:]]
    def _set_xyz_cartn(self, xyzList):
        for atom,xyz_cartn in zip(self, xyzList):
            atom.xyz_cartn = xyz_cartn
    xyz_cartn = property(_get_xyz_cartn, _set_xyz_cartn, doc =
        """absolute Cartesian coordinates of all atoms""" )   
    
    # symbols
    def _get_symbols(self):
        return [atom.symbol for atom in self[:]]
    def _set_symbols(self, symbolList):
        for atom,symbol in zip(self, symbolList):
            atom.symbol = symbol
    symbols = property(_get_symbols, _set_symbols, doc =
        """symbols of all atoms""" )  
    
    # forces
    def _get_forces(self):
        return [atom.force for atom in self]
    def _set_forces(self, forceList):
        for atom,force in zip(self, forceList):
            atom.force = force
    forces = property(_get_forces, _set_forces, doc =
        """forces on all atoms""" )   
    
    # charges
    def _get_charges(self):
        return [atom.charge for atom in self]
    def _set_charges(self, chargeList):
        for atom,charge in zip(self, chargeList):
            atom.charge = charge
    charges = property(_get_charges, _set_charges, doc =
        """charges on all atoms in electron units""" )   
    
################################################################################################    
# geometry and symmetry methods--these should be farmed out to Geometry class which does all this--see vimm
################################################################################################ 

    def distance(self, id0, id1):
        """Distance between 2 atoms, no periodic boundary conditions.

        id0 -- zero based index of the first atom or a string label
               such as "Na1"
        id1 -- zero based index or string label of the second atom.

        Return float.
        Raise ValueError for invalid arguments.
        """
        a0 = self.getAtom(id0)
        a1 = self.getAtom(id1)
        return self.lattice.dist(a0.xyz, a1.xyz)


    def angle(self, a0, a1, a2):
        """angle at atom a1 in degrees"""
        u10 = a0.xyz - a1.xyz
        u12 = a2.xyz - a1.xyz
        return self.lattice.angle(u10, u12)


    def symConsistent(self, decimal=7):
        from numpy.testing import assert_array_almost_equal
        verdict = True
        nonCompliantAtomPosition = None
        nonCompliantSymOp = None
        def arrayInList(trialArray,arrayList):
            matchesQ=False
            for certainArray in arrayList:
                if (numpy.round(trialArray-certainArray, decimal)==0).all():
                    matchesQ = True
            return matchesQ
        sg = self.sg
        # as long as the structure is given a lattice, each atom.xyz should be in fractional coords
        atomPositions = [atom.xyz for atom in self[:]]
        for atomPosition in atomPositions:
            for symop in sg.symop_list:
                rawPosition = numpy.dot(symop.R, atomPosition) + symop.t
                fracPos, intPos = numpy.modf(numpy.round(rawPosition, decimal))
                newPosition = numpy.mod(fracPos,1)
                if not arrayInList(newPosition, atomPositions):
                    verdict = False
                    nonCompliantAtomPosition = atomPosition
                    nonCompliantSymOp = symop
                    break
        return verdict, nonCompliantAtomPosition, nonCompliantSymOp
    
    def placeInLattice(self, new_lattice):
        """place structure into new_lattice coordinate system

        sets lattice to new_lattice and recalculate fractional coordinates
        of all atoms so their absolute positions remain the same

        return self
        """
        Tx = numpy.dot(self.lattice.base, new_lattice.recbase)
        Tu = numpy.dot(self.lattice.normbase, new_lattice.recnormbase)
        for a in self:
            a.xyz = list(numpy.dot(a.xyz, Tx))
        self.lattice = new_lattice
        return self
    
    def computeDistances(self, maxdist=30, latticeRange=[2,2,2]):
        """ unitcell.computeDistances(self, [nx,ny,nz]):
        builds up a Big multiple dictionary, namely
        self.distances[atom1][atom2][(DX,DY,DZ)]
        (DX,DY,DZ) are integer numbers specifying the cell containing atom2,
        relatively to atom1.
        DX,DY,DZ run from -nx to nx, -ny to ny, -nz to nz, respectively."""
        distances = {}
        idlist = self.getLabels()
        #idlist = self.getIds()
        for iA in range(0, len(idlist)):
            idA = idlist[iA]
            distances[idA] = {}
            for iB in range(0, len(idlist)):
                idB = idlist[iB]
                distances[idA][idB]={}
                for tx in range(-latticeRange[0],latticeRange[0]+1):
                    for ty in range(-latticeRange[1],latticeRange[1]+1):
                        for tz in range(-latticeRange[2],latticeRange[2]+1):
                            posA = self.getCartesianPosition(idA)
                            posB = self.getCartesianPosition(idB) + numpy.dot([tx,ty,tz], self._lattice)
                            dist = numpy.sqrt(numpy.sum( (posB-posA) * (posB-posA) ))
                            if(dist<maxdist):
                                distances[idA][idB][(tx,ty,tz)] = dist
        self._distances = distances
        return distances

##     def findTetrahedra(self, list4ids, latticeVectors=[(0,0,0),(0,0,0),(0,0,0),(0,0,0)]):
##         """Searches the tetrahedra that equivalent by symmetry to the reference tetrahedron
##         defined by the list4ids (4 Ids) and the cells indices where vertices lie."""        
##         # from OpenPhonon:
##         # OP_cella. FindTetragons(self, Latoms, CellsCoo=[(0,0,0),(0,0,0),(0,0,0),(0,0,0)])
##         # This is the key routine to find good candidates to symmetry group.
##         # Latoms is a list of the kind [  (Aa,kA) ,(Bb,kB)....] 
##         # ( i.e. couples formed by atom-name and position in the position list.
##         # Latoms must be formed of four  atoms defining a non-degenerate tetraedron.
##         # A check is permormed on the non-degeneracy.
##         # The funtions finds all the possible equivalent tetraedrons (which have the same
##         # set of distances, and the same atom kinds)
##         # The function return the list of all these tetraedrons
        
##         if type(list4ids) is not type([]):
##             raise ValueError, 'list4ids should be a list of 4 site Ids.'
##         if len(list4ids) != len(latticeVectors):
##             raise ValueError, 'There should be as many site Ids as lattice vectors.'
##         if len(list4ids) != 4:
##             raise ValueError, 'Need 4 sites to define a tetrahedron.'
        
##         tetraVertices = [self.cartesianPositionInLattice(id,lattvec)
##                          for (id,lattvec) in zip(list4ids,latticeVectors)]
##         # compute vectors for edges of tetrahedron from first point:
##         edgeVectors = tetraVertices[1:4] - tetraVertices[0]

##         # check for non-degeneracy:
##         det = la.det(edgeVectors)
##         if(abs(det) < 1e-6):
##             raise ValueError, 'determinant smaller than 1e-6: degenerate reference tetrahedron.'

###########################
# IO
###########################
    def read(self, filename, format='auto'):
        """Load structure from a file, any original data may become lost.

        filename -- file to be loaded
        format   -- all structure formats are defined in Parsers submodule,
                    when format == 'auto' all Parsers are tried one by one

        Return instance of data Parser used to process file.  This
        can be inspected for information related to particular format.
        """
        from .Parsers import getParser
        p = getParser(format)
        new_structure = p.parseFile(filename)
        # reinitialize data after successful parsing
        # avoid calling __init__ from a derived class
        Structure.__init__(self)
        if new_structure is not None:
            self.__dict__.update(new_structure.__dict__)
            self[:] = new_structure
        if not self.description:
            self.generateDescription()
#            import os.path
#            tailname = os.path.basename(filename)
#            tailbase = os.path.splitext(tailname)[0]
#            self.description = tailbase
        return p

    def readStr(self, s, format='auto'):
        """Load structure from a string, any original data become lost.

        s        -- string with structure definition
        format   -- all structure formats are defined in Parsers submodule,
                    when format == 'auto' all Parsers are tried one by one

        Return instance of data Parser used to process input string.  This
        can be inspected for information related to particular format.
        """
        from .Parsers import getParser
        p = getParser(format)
        new_structure = p.parse(s)
        # reinitialize data after successful parsing
        # avoid calling __init__ from a derived class
        Structure.__init__(self)
        if new_structure is not None:
            self.__dict__.update(new_structure.__dict__)
            self[:] = new_structure
        return p

    def write(self, filename, format, **kwds):
        """Save structure to file in the specified format

        No return value.

        Note: available structure formats can be obtained by:
            from Parsers import formats
        """
        from .Parsers import getParser
        p = getParser(format)
        p.filename = filename
        s = p.tostring(self, **kwds)
        f = open(filename, 'wb')
        f.write(s)
        f.close()
        return

    def writeStr(self, format, **kwds):
        """return string representation of the structure in specified format

        Note: available structure formats can be obtained by:
            from Parsers import formats
        """
        from .Parsers import getParser
        p = getParser(format)
        s = p.tostring(self, **kwds)
        return s

    ##############################################################################
    # overloaded list methods
    ##############################################################################

    def append(self, a, copy=True):
        """Append atom to a structure and update its lattice attribute.

        a    -- instance of Atom
        copy -- flag for appending a copy of a.
                When False, append a and update a.owner.

        No return value.
        """
        self._uncache('labels')
        adup = copy and Atom(a) or a
        adup.lattice = self.lattice
        list.append(self, adup)
        return

    def insert(self, idx, a, copy=True):
        """Insert atom a before position idx in this Structure.

        idx  -- position in atom list
        a    -- instance of Atom
        copy -- flag for inserting a copy of a.
                When False, append a and update a.lattice.

        No return value.
        """
        self._uncache('labels')
        adup = copy and Atom(a) or a
        adup.lattice = self.lattice
        list.insert(self, idx, adup)
        return

    def extend(self, atoms, copy=True):
        """Extend Structure by appending copies from a list of atoms.

        atoms -- list of Atom instances
        copy  -- flag for extending with copies of Atom instances.
                 When False extend with atoms and update their lattice
                 attributes.

        No return value.
        """
        self._uncache('labels')
        if copy:    adups = [Atom(a) for a in atoms]
        else:       adups = atoms
        for a in adups: a.lattice = self.lattice
        list.extend(self, adups)
        return

    def __setitem__(self, idx, a, copy=True):
        """Set idx-th atom to a.

        idx  -- index of atom in this Structure
        a    -- instance of Atom
        copy -- flag for setting to a copy of a.
                When False, set to a and update a.lattice.

        No return value.
        """
        self._uncache('labels')
        adup = copy and Atom(a) or a
        adup.lattice = self.lattice
        list.__setitem__(self, idx, adup)
        return

    def __setslice__(self, lo, hi, atoms, copy=False):
        """Set Structure slice from lo to hi-1 to the sequence of atoms.

        lo    -- low index for the slice
        hi    -- high index of the slice
        atoms -- sequence of Atom instances
        copy  -- flag for using copies of Atom instances.  When False, set
                 to existing instances and update their lattice attributes.

        No return value.
        """
        self._uncache('labels')
        if copy:    
            adups = [Atom(a) for a in atoms]
        else:       
            adups = atoms
        for a in adups: a.lattice = self.lattice
        list.__setslice__(self, lo, hi, adups)

    ####################################################################
    # property handlers
    ####################################################################

    # lattice

    def _get_lattice(self):
        if not hasattr(self, '_lattice'):
            self._lattice = Lattice()
        return self._lattice
    def _set_lattice(self, value):
        for a in self:  a.lattice = value
        self._lattice = value
        return
    lattice = property(_get_lattice, _set_lattice, doc =
        "Coordinate system for this Structure.")
    
    # space group

    def _get_spaceGroup(self):
        if not self._sg:
            from .SpaceGroups import GetSpaceGroup
            self._sg = GetSpaceGroup(self._sgid)
        return self._sg

    def _set_spaceGroup(self, item):
        from .SpaceGroups import SpaceGroup
        if isinstance(item, SpaceGroup):
            self._sg = item
            self._sgid = item.number
        else:
            from .SpaceGroups import GetSpaceGroup
            self._sg = GetSpaceGroup(item)
            self._sgid = item

    sg = property(_get_spaceGroup, _set_spaceGroup, doc =
        """Space group for this structure.  This can be set 
        by instantiating with a new spacegroup class or with a space group id.
        One can also use the explicit setter.""")
    
    ####################################################################
    # additional setter
    ####################################################################

    def setSg(self, sgid):
        self.sg = sgid
        return

    ####################################################################
    # protected methods
    ####################################################################

    def _update_labels(self):
        """Update the _labels dictionary of unique string labels of atoms.

        No return value.
        """
        kv = list(zip(self.getLabels(), self[:]))
        self._labels = dict(kv)
        self._labels_cached = True
        return


    def _uncache(self, *args):
        """Reset cached flag for a list of internal attributes.

        *args -- list of strings, currently supported are "labels"

        No return value.
        Raise AttributeError for any invalid args.
        """
        for a in args:
            attrname = "_" + a + "_cached"
            setattr(self, attrname, False)
        return



