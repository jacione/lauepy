__doc__ = """Test conversion of CIF file into a UnitCell instance."""

from Translator import Translator
from UnitCell import UnitCell
from Atom import Atom
import math
from math import cos,sin,pi

fileInString = "corundum.cif"

def MInv(uc):
    """gives the Cartesian normalization with respect the x axis of an arbitrary unit cell uc input as a tuple. Implicit assumption is that a_vector is along x and b_vector in the (x,y) plane."""
    conv=math.pi/180.
    # a=uc[0];b=uc[1];c=uc[2];al=uc[3];be=uc[4];ga=uc[5]
    (a,b,c,al,be,ga)=uc
    cosAlStar=(math.cos(conv*be)*math.cos(conv*ga) - math.cos(conv*al))/(math.sin(conv*be)*math.sin(conv*ga))
    V = a*b*c*(1 - math.cos(conv*al)**2. - math.cos(conv*be)**2. - math.cos(conv*ga)**2. + 2.*math.cos(conv*al)*math.cos(conv*be)*math.cos(conv*ga))**(1/2.)
    cStar = a*b*math.sin(conv*ga)/V
    return [[a, 0, 0], [b*math.cos(conv*ga), b*math.sin(conv*ga), 0], [c*math.cos(conv*be), -c*math.sin(conv*be)*cosAlStar, 1/cStar]]


# Brandon: I couldn't find out to retrieve the unit cell vectors
# from the CIF file using your Translator...
# This is a temporary fix.
# Olivier
from pyparsing import *

number = Combine( Optional('-') + ( '0' | Word('123456789',nums) ) + \
                  Optional( '.' + Word(nums) ) + \
                  Optional( Word('eE',exact=1) + Word(nums+'+-',nums) ) )

def convertNumbers(s,l,toks):
    n = toks[0]
    try: return int(n)
    except ValueError as ve: return float(n)
    raise

number.setParseAction( convertNumbers )

cell_length_a = (Literal("_cell_length_a").setParseAction( replaceWith("a")) + number)
cell_length_b = (Literal("_cell_length_b").setParseAction( replaceWith("b")) + number)
cell_length_c = (Literal("_cell_length_c").setParseAction( replaceWith("c")) + number)

cell_angle_alpha = (Literal("_cell_angle_alpha").setParseAction( replaceWith("alpha")) + number)
cell_angle_beta = (Literal("_cell_angle_beta").setParseAction( replaceWith("beta")) + number)
cell_angle_gamma = (Literal("_cell_angle_gamma").setParseAction( replaceWith("gamma")) + number)

cell_parsers = [cell_length_a, cell_length_b, cell_length_c, cell_angle_alpha, cell_angle_beta, cell_angle_gamma]

try:
    ifile = open(fileInString, 'r')
except IOError as xxx_todo_changeme:
    (errno, strerror) = xxx_todo_changeme.args
    print(("I/O error(%s): %s" % (errno, strerror)))

string = ifile.readlines()

cell_params = []
cell_dict = {}
for parser in cell_parsers:
    for data, dataStart, dataEnd in parser.scanString(string):
        cell_params.append(data[1])
        cell_dict[data[0]]=data[1]

ifile.close()

### end of fix

t=Translator()
t.filenameIn = fileInString
ciflist = t.cifToAtomAndCoordinateList()

cellvectors = MInv(cell_params)

uc = UnitCell(cellvectors)

for item in ciflist:
    X = item[0]
    pos = [float(x) for x in item[1:4]]
    uc.addAtom(Atom(symbol=X), pos, "")
    # Note: this only works if elements' symbols in CIF are "Ab" format

