# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 16:27:31 2016

@author: ruqingxu
"""

import scipy.constants as constants

def keV_to_inv_m(x):
    return x*1000*constants.e/constants.hbar/constants.c

def inv_m_to_keV(x):
    return x*constants.c*constants.hbar/constants.e/1000
