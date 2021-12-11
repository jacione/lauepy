# -*- coding: utf-8 -*-
"""
Collection of classes & functions for crystal-lattice-related handling

External dependency as of Jan 2016:
    "matter" package from https://github.com/danse-inelastic/matter
                 also see http://dev.danse.us/trac/inelastic/wiki/crystal

To do as of Jan 2016:
    * Xtal class
    * enable x-ray structure factor evaluation from reference data

Xu, Ruqing
Created Jan 2016
"""

import numpy as np
from .extutils import matter as __matter
import periodictable as PT

mAtom = __matter.Atom
mLattice = __matter.Lattice

class AtomInCell(object):
    '''
    a simpler class that wraps over matter.Atom;
    defines atom type & fractional position
    default is a Si atom at origin
    '''
    def __init__(self,symbol='Si',xf=.0,yf=.0,zf=.0):
        self.symbol = symbol
        self.xyzf = [xf,yf,zf]
        self._matom_obj = mAtom(symbol,[xf,yf,zf])
    @property
    def xf(self):
        return self.xyzf[0]
    @property
    def yf(self):
        return self.xyzf[1]
    @property
    def zf(self):
        return self.xyzf[2]
    @property
    def Z(self):
        return self._matom_obj.Z
    def __repr__(self):
        return str((self.symbol,self.xyzf))
    def __str__(self):
        return str((self.symbol,self.xyzf))

class Xtal(object):
    '''
    Note: All SI units
    
    todo: 
        * documentation
        * .loadfromXML()
        * .get_structure_factor()
    '''
    def __init__(self,a,b,c,alpha,beta,gamma,atomlist=[],pointgroup=0):
        self.abc = [a,b,c]
        self.angles = [alpha,beta,gamma]
        if atomlist:
            self.atoms = atomlist
        else:   # default is a single Si atom
            self.atoms = [AtomInCell()] 
        self.pointgroup = pointgroup # for future use
        self.description = '' # one can add description here
        self._mlatt_obj = mLattice(a*1e10,b*1e10,c*1e10,alpha,beta,gamma)
        self._mlatt_rec = self._mlatt_obj.reciprocal()
        
    @property
    def num_atoms(self):
        return len(self.atoms)
    
    @property
    def real_bases(self):
        '''
        Gives 3 real lattice basis vectors, in SI units.
        Stored as row vectors, e.g., result[0] is vector a
        '''
        vec = self._mlatt_obj.base
        return np.array(vec) * 1.0e-10
    
    @property
    def rec_bases(self):
        '''
        Gives 3 reciprocal lattice basis vectors, in SI units.
        Stored as row vectors, e.g., result[0] is vector a_star
        '''
        return np.array(self._mlatt_rec.base) * 2*np.pi * 1.0e10
    
    @property
    def abc_stars(self):
        params = self._mlatt_rec.abcABG()
        return np.array(params[0:3]) * np.pi * 2.0e10
    @property
    def angles_rec(self):
        params = self._mlatt_rec.abcABG()
        return np.array(params[3:6])
    
    def get_coord(self,ijk):
        '''
        generate coordinates of a given ijk vector
        '''
        return np.array(self._mlatt_obj.cartesian(ijk)) * 1.0e-10
    
    def get_q_vec(self,hkl):
        '''
        generate coordinates of q vector from a given hkl vector
        '''
        return np.array(self._mlatt_rec.cartesian(hkl)) * 2*np.pi * 1.0e10    
    
    def get_q_mag(self,hkl):
        '''
        compute magnitude of q vector from a given hkl vector
        '''
        q_vec = self.get_q_vec(hkl) 
        return np.linalg.norm(q_vec)

    def get_atom_coord(self,i=-1):
        '''
        i is atom index, if not provided, return coordinates of all,
        in which case, result is Nx3 matrix; otherwise is a 3-vector
        '''
        if i >= 0:
            return self.get_coord(self.atoms[i].xyzf)
        else:
            result = [self.get_coord(atom.xyzf) for atom in self.atoms]
            return np.array(result)
        
    def get_structure_factor(self,hkl):
        '''
        currently use Z of each atom,
        
        todo: use x-ray reference data for more real calculation
        '''
        q = self.get_q_vec(hkl)
        q_len_ang = np.linalg.norm(q) * 1.0e-10  # length of q vector in inverse angstrom
        sf = 0.0j
        for atom in self.atoms:
            ## compute atomic scattering factor ##
            element = PT.elements[atom.Z]
            asf = element.xray.f0(q_len_ang)
            if np.isnan(asf):
                asf = 0.0
            ## assemble structure factor ##
            atom_pos = self.get_coord(atom.xyzf)
            phase = np.dot(atom_pos,q)
            sf += np.exp(1j*phase)*asf

        return sf.real**2 + sf.imag**2
    
    def reflection_is_allowed(self,hkl):
        '''
        Check whether a particular reflection (h,k,l) is allowed
        in reciprocal space.
        '''
        sf = self.get_structure_factor(hkl)
        return np.absolute(sf) > 1.0e-6
        
    def loadfromfile(self):
        pass