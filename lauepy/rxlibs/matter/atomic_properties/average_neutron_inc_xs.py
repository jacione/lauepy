
__doc__ = """Average neutron incoherent scattering cross-section for the elements
(average over isotopes in their natural abundancies) in barns,
as provided on the NIST website."""

def lookup( atom ):
    return avg_neutron_inc_xs[ atom.Z ]

avg_neutron_inc_xs = [
    None, # X
    80.2600, # H
    0.0000, # He
    0.9200, # Li
    0.0018, # Be
    1.7000, # B
    0.0010, # C
    0.5000, # N
    0.0008, # O
    0.0008, # F
    0.0080, # Ne
    1.6200, # Na
    0.0800, # Mg
    0.0082, # Al
    0.0040, # Si
    0.0050, # P
    0.0070, # S
    5.3000, # Cl
    0.2250, # Ar
    0.2700, # K
    0.0500, # Ca
    4.5000, # Sc
    2.8700, # Ti
    5.0800, # V
    1.8300, # Cr
    0.4000, # Mn
    0.4000, # Fe
    4.8000, # Co
    5.2000, # Ni
    0.5500, # Cu
    0.0770, # Zn
    0.1600, # Ga
    0.1800, # Ge
    0.0600, # As
    0.3200, # Se
    0.1000, # Br
    0.0100, # Kr
    0.5000, # Rb
    0.0600, # Sr
    0.1500, # Y
    0.0200, # Zr
    0.0024, # Nb
    0.0400, # Mo
    0.5000, # Tc
    0.4000, # Ru
    0.3000, # Rh
    0.0930, # Pd
    0.5800, # Ag
    3.4600, # Cd
    0.5400, # In
    0.0220, # Sn
    0.0070, # Sb
    0.0900, # Te
    0.3100, # I
    0.0000, # Xe
    0.2100, # Cs
    0.1500, # Ba
    1.1300, # La
    0.0010, # Ce
    0.0150, # Pr
    9.2000, # Nd
    1.3000, # Pm
    39.0000, # Sm
    2.5000, # Eu
    151.0000, # Gd
    0.0040, # Tb
    54.4000, # Dy
    0.3600, # Ho
    1.1000, # Er
    0.1000, # Tm
    4.0000, # Yb
    0.7000, # Lu
    2.6000, # Hf
    0.0100, # Ta
    1.6300, # W
    0.9000, # Re
    0.3000, # Os
    0.0000, # Ir
    0.1300, # Pt
    0.4300, # Au
    6.6000, # Hg
    0.2100, # Tl
    0.0030, # Pb
    0.0084, # Bi
    None , # Po
    None , # At
    None , # Rn
    None , # Fr
    0.0000, # Ra
    None , # Ac
    0.0000, # Th
    0.1000, # Pa
    0.0050, # U
    0.5000, # Np
    None , # Pu
    0.3000, # Am
    None , # Cm
    None , # Bk
    None , # Cf
    None , # Es
    None , # Fm
    None , # Md
    None , # No
    None] # Lw

# end of averaged_neutron_incoh_xs.py

