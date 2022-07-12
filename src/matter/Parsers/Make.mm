# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2005  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJECT = matter
PACKAGE = Parsers

#--------------------------------------------------------------------------
#

BUILD_DIRS = \
	
OTHER_DIRS = \

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)

#--------------------------------------------------------------------------
#

all: export
	BLD_ACTION="all" $(MM) recurse


#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
	P_auto.py \
	P_bratoms.py \
	P_cif.py \
	P_pdb.py \
	P_pdffit.py \
	P_rawxyz.py \
	P_xcfg.py \
	P_xyz.py \
	__init__.py \
	parser_index_mod.py \
	structureparser.py \


export:: export-package-python-modules

# version
# $Id$

# End of file
