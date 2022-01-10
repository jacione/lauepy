
__doc__ = """Average neutron incoherent scattering cross-section for the elements
(average over isotopes in their natural abundancies) in barns,
as provided on the NIST website."""

def lookup( atom ):
    return avg_neutron_abs_xs[ atom.Z ]

avg_neutron_abs_xs = [
    None, # X
    0.33260, # H
    0.00747, # He
    70.50000, # Li
    0.00760, # Be
    767.00000, # B
    0.00350, # C
    1.90000, # N
    0.00019, # O
    0.00960, # F
    0.03900, # Ne
    0.53000, # Na
    0.06300, # Mg
    0.23100, # Al
    0.17100, # Si
    0.17200, # P
    0.53000, # S
    33.50000, # Cl
    0.67500, # Ar
    2.10000, # K
    0.43000, # Ca
    27.50000, # Sc
    6.09000, # Ti
    5.08000, # V
    3.05000, # Cr
    13.30000, # Mn
    2.56000, # Fe
    37.18000, # Co
    4.49000, # Ni
    3.78000, # Cu
    1.11000, # Zn
    2.75000, # Ga
    2.20000, # Ge
    4.50000, # As
    11.70000, # Se
    6.90000, # Br
    25.00000, # Kr
    0.38000, # Rb
    1.28000, # Sr
    1.28000, # Y
    0.18500, # Zr
    1.15000, # Nb
    2.48000, # Mo
    20.00000, # Tc
    2.56000, # Ru
    144.80000, # Rh
    6.90000, # Pd
    63.30000, # Ag
    2520.00000, # Cd
    193.80000, # In
    0.62600, # Sn
    4.91000, # Sb
    4.70000, # Te
    6.15000, # I
    23.90000, # Xe
    29.00000, # Cs
    1.10000, # Ba
    8.97000, # La
    0.63000, # Ce
    11.50000, # Pr
    50.50000, # Nd
    168.40000, # Pm
    5922.00000, # Sm
    4530.00000, # Eu
    49700.00000, # Gd
    23.40000, # Tb
    994.00000, # Dy
    64.70000, # Ho
    159.00000, # Er
    100.00000, # Tm
    34.80000, # Yb
    74.00000, # Lu
    104.10000, # Hf
    20.60000, # Ta
    18.30000, # W
    89.70000, # Re
    16.00000, # Os
    425.00000, # Ir
    10.30000, # Pt
    98.65000, # Au
    372.30000, # Hg
    3.43000, # Tl
    0.17100, # Pb
    0.03380, # Bi
    None, # Po
    None, # At
    None, # Rn
    None, # Fr
    12.80000, # Ra
    None, # Ac
    7.37000, # Th
    200.60000, # Pa
    7.57000, # U
    175.90000, # Np
    None, # Pu
    75.30000, # Am
    None, # Cm
    None, # Bk
    None, # Cf
    None, # Es
    None, # Fm
    None, # Md
    None, # No
    None] # Lw
    
# end of averaged_neutron_abs_xs.py

