##############################################################################
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Index of recognized structure formats, their IO capabilities and
associated modules where they are defined.  Plugins for new structure
formats need to be added to the parser_index dictionary in this module.
"""

__id__ = "$Id: parser_index_mod.py 3032 2009-04-08 19:15:37Z juhas $"

parser_index = {

    # automatic format detection - tries all parsers one by one
    'auto' : {
        'module' : 'P_auto',
        'file_extension' : '',
        'file_pattern' : '*.*',
        'has_input' : True,
        'has_output' : False,
        },

    # CIF format
    'cif' : {
        'module' : 'P_cif',
        'file_extension' : '.cif',
        'file_pattern' : '*.cif',
        'has_input' : True,
        'has_output' : True,
        },

    # PDB format
    'pdb' : {
        'module' : 'P_pdb',
        'file_extension' : '.pdb',
        'file_pattern' : '*.pdb',
        'has_input' : True,
        'has_output' : True,
        },

    # standard xyz file
    'xyz' : {
        'module' : 'P_xyz',
        'file_extension' : '.xyz',
        'file_pattern' : '*.xyz',
        'has_input' : True,
        'has_output' : True,
        },

    # raw xyz file (element labels optional)
    'rawxyz' : {
        'module' : 'P_rawxyz',
        'file_extension' : '.xyz',
        'file_pattern' : '*.xyz',
        'has_input' : True,
        'has_output' : True,
        },

    # AtomEye extended configuration format
    'xcfg' : {
        'module' : 'P_xcfg',
        'file_extension' : '',
        'file_pattern' : '*.xcfg|*.eye|*.cfg',
        'has_input' : True,
        'has_output' : True,
        },
        
        
    # PDFfit structure format
    'pdffit' : {
        'module' : 'P_pdffit',
        'file_extension' : '.stru',
        'file_pattern' : '*.stru|*.rstr',
        'has_input' : True,
        'has_output' : True,
        },
        
        # Bruce Ravel's atoms format
#    'bratoms' : {
#        'module' : 'P_bratoms',
#        'file_extension' : '.inp',
#        'file_pattern' : '*.inp',
#        'has_input' : True,
#        'has_output' : True,
#        },
}

# End of file
