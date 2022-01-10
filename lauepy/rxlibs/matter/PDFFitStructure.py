##############################################################################
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""definition of PDFFitStructure class derived from Structure
"""

__id__ = "$Id: pdffitstructure.py 3458 2009-07-08 06:38:31Z juhas $"

from danse.ins.matter import Structure

##############################################################################
class PDFFitStructure(Structure):
    """PDFFitStructure --> Structure with extra pdffit member

    Data members:
        pdffit -- dictionary for storing following extra parameters from
                  PDFFit structure files:
                      'scale', 'delta1', 'delta2', 'sratio',
                      'rcut', 'spcgr', 'dcell', 'ncell'
    """

    def __init__(self, *args, **kwargs):
        Structure.__init__(self, *args, **kwargs)
        # do not overwrite self.pdffit, when instantiated as a copy
        # of existing PDFFitStructure
        if not hasattr(self, 'pdffit'):
            self.pdffit = {
                'scale' : 1.0,
                'delta1' : 0.0,
                'delta2' : 0.0,
                'sratio' : 1.0,
                'rcut' : 0.0,
                'spcgr' : 'P1',
                'spdiameter' : 0.0,
                'stepcut' : 0.0,
                'dcell' : 6*[0.0],
                'ncell' : [1, 1, 1, 0],
            }
        return

    def read(self, filename, format='auto'):
        """Same as Structure.read, but update spcgr value in
        self.pdffit when parser can get spacegroup.

        Return instance of StructureParser used to load the data.
        See Structure.read() for more info.
        """
        p = Structure.read(self, filename, format)
        sg = getattr(p, 'spacegroup', None)
        if sg:  self.pdffit['spcgr'] = sg.short_name
        return p

    def readStr(self, s, format='auto'):
        """Same as Structure.readStr, but update spcgr value in
        self.pdffit when parser can get spacegroup.

        Return instance of StructureParser used to load the data.
        See Structure.readStr() for more info.
        """
        p = Structure.readStr(self, s, format)
        sg = getattr(p, 'spacegroup', None)
        if sg:  self.pdffit['spcgr'] = sg.short_name
        return p

# End of PDFFitStructure
