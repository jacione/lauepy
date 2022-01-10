##############################################################################
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Exceptions used in Structure package.
"""

__id__ = "$Id: StructureErrors.py 2825 2009-03-09 04:33:12Z juhas $"

class StructureFormatError(Exception):
    """Exception for failed IO from Structure file
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

# End of StructureFormatError


class LatticeError(Exception):
    """Exception for impossible lattice parameters.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

# End of LatticeError


class SymmetryError(Exception):
    """Exception raised for invalid symmetry operations.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

# End of SymmetryError


class IsotropyError(Exception):
    """Exception raised for invalid operations on isotropic atoms.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

# End of IsotropyError
