# This algorithm was taken from Campos/ASE
# and updated to use numpy instead of numeric.
# Olivier Delaire (4/2007)

__doc__ = """Generates a Monkhorst-Pack grid of points of order (n1,n2,n3).
The points are generated in the open cube
]-1/2,1/2[ x ]-1/2,1/2[ x ]-1/2,1/2[/"""

import numpy as np

def MonkhorstPack(size):
    """Generates a Monkhorst-Pack grid of points of order (n1,n2,n3).
    The points are generated in the open cube
    ]-1/2,1/2[ x ]-1/2,1/2[ x ]-1/2,1/2[/"""
    
    if (np.size(size) != 3):
        raise ValueError("Monkhorst-Pack grid size argument must have 3 values for 3D grid.")
    elif (0 in size):
        raise ValueError("Monkhorst-Pack grid order along any dimension cannot be zero.")        
    else:
        kpts = np.swapaxes(np.indices(size, np.float), 0, 3)
        kpts = np.reshape(kpts, (-1, 3))
        return (kpts + (0.5, 0.5, 0.5)) / size - (0.5, 0.5, 0.5)

if __name__ == '__main__':
    print((MonkhorstPack((1, 1, 1))))
