#__import__('pkg_resources').declare_namespace(__name__)

import sys
vinfo =  sys.version_info
if vinfo[0] != 3: raise NotImplementedError
METACLASS_CTOR_NEEDS_TYPE = vinfo[1]>5


from . import crystalIO
from . import atomic_properties
from .StructureErrors import StructureFormatError,LatticeError,SymmetryError,IsotropyError
from .Atom import Atom
from .Lattice import Lattice
from .Structure import Structure

#def atom(**kwds):
#    from .Atom import Atom
#    return Atom(**kwds)
#
#def lattice(**kwds):
#    from .Lattice import Lattice
#    return Lattice(**kwds)
#
#def structure(**kwds):
#    from .Structure import Structure
#    return Structure(**kwds)
    
    
    



# obtain version information
from .version import __version__


