__doc__ = """small methods for 3d plotting"""

def plotUnitCell(unitcell):
    """Plots a unit cell in 3D."""
    from ExcitationSlicer.vtktools.vtkplotting import VTKPlotAtoms
    from crystal.crystalIO.converters import unitCell2ListOfAtom
    loa = unitCell2ListOfAtom(unitcell)
    plot = VTKPlotAtoms(loa)
    #raw_input("Press Enter.")    
    return plot
