
TriposMOL2AtomTypes = {
#Carbon (C):

    'C.3'  : 'C:Tet' ,   #.Carbon in sp³ hybridization (tetrahedral geometry, single bonds, such as in alkanes).
    'C.2'  : 'C:Tri' ,   #.Carbon in sp² hybridization (planar geometry, double bonds, such as in alkenes).
    'C.1'  : 'C:Lin' ,   #.Carbon in sp hybridization (linear geometry, triple bonds, such as in alkynes).
    'C.ar' : 'C:Res' ,   #.Carbon in aromatic systems (part of a conjugated pi system, such as in benzene).
    'C.cat': 'C:Tri' ,   #.Carbon in a cationic state (carbocation, typically sp² hybridized).

#Nitrogen (N):

    'N.3':   'N:Tet' ,   #.Nitrogen in sp³ hybridization (tetrahedral geometry, single bonds, as in amines).
    'N.2':   'N:Tri' ,   #.Nitrogen in sp² hybridization (planar geometry, double bonds, such as in imines).
    'N.1':   'N:Lin' ,   #.Nitrogen in sp hybridization (linear geometry, triple bonds, such as in nitriles).
    'N.am':  'N:Res' ,   #.Nitrogen in an amide group (sp² hybridized, as in peptide bonds).
    'N.pl3': 'N:Tri' ,   #.Nitrogen in a planar, sp² hybridized state, not carrying a formal charge (such as in aromatic heterocycles like pyridine).
    'N.4':   'N:Tet' ,   #.Nitrogen in a sp³ hybridized state, carrying a formal positive charge (such as in quaternary ammonium groups).

#Oxygen (O):

    'O.3'  : None ,    #.Oxygen in sp³ hybridization (tetrahedral geometry, single bonds, such as in alcohols and ethers).
    'O.2'  : None ,    #.Oxygen in sp² hybridization (planar geometry, double bonds, such as in carbonyl groups, C=O).
    'O.co2': None ,    #.Oxygen in a carboxylate group (charged, as in carboxylate anions, COO⁻).

#Sulfur (S):

    'S.3' :  None ,    #.Sulfur in sp³ hybridization (tetrahedral geometry, single bonds, such as in sulfides).
    'S.2' :  None ,    #.Sulfur in sp² hybridization (planar geometry, such as in thioethers and thioketones).
    'S.o' :  None ,    #.Sulfur in a sulfoxide group (S=O).
    'S.o2':  None ,    #.Sulfur in a sulfone group (O=S=O).
    'S.3+':  None ,    #.Sulfur in a positively charged sp³ hybridized state (such as in sulfonium ions).

#Phosphorus (P):

    'P.3' : None  ,    #.Phosphorus in sp³ hybridization (tetrahedral geometry, such as in phosphines).
    'P.3+': None  ,    #.Phosphorus in a positively charged sp³ hybridized state (such as in phosphonium ions).

#Halogens:

    'F' : None ,  #.Fluorine.
    'Cl': None ,  #.Chlorine.
    'Br': None ,  #.Bromine.
    'I' : None ,  #.Iodine.

#Metals:

    'Li': None, #Lithium.
    'Na': None, #Sodium.
    'K' : None, #Potassium.
    'Ca': None, #Calcium.
    'Mg': None, #Magnesium.
    'Zn': None, #Zinc.
    'Fe': None, #Iron.
    'Cu': None, #Copper.
    'Mn': None, #Manganese.

#Hydrogen (H):

    'H': None, #Hydrogen atom.
}
'''
Special Atom Types:

    Du: Dummy atom (used in special cases for force field calculations or placeholders).
    LP: Lone pair of electrons (not typically used in many force fields).
    Any: Wildcard atom type (used in SMARTS and substructure searching).
'''
