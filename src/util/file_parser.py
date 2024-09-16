from pprint import pprint


def get_file_type(filename):
    file_type = filename.split('.')
    return file_type[-1]
    

def read_SIMPLE_txt (filein = None):
    """ Function doc """    
    charges = []
    
    with open(filein, 'r') as _file:
        filetext = _file.readlines()
        
        for line in filetext:
            print(line)
            #line = line.strip()
            if len(line) > 0:
                try:
                    line = float(line)
                    charges.append(line)
                except:
                    charges = False
                    #break
            else:
                pass
    return charges
    
    
def read_MOPAC_out (filein = None):
    """ Function doc """

def read_MOPAC_aux (filein = None):
    """ Function doc """
    data = open(filein, 'r')
    data = data.read()
    
    #.ATOM_CORE
    try:
        data2 = data.split('ATOM_CORE')
        data2 = data2[-1].split('ATOM_X:ANGSTROMS')
        data2 = data2[0] 

        data2 = data2.split('=')
        data2 = data2[-1].split()

        ATOM_CORE = data2
    except:
        ATOM_CORE = None

    #.HEAT_OF_FORMATION
    try:
        data2 = data.split('HEAT_OF_FORMATION:KCAL/MOL')
        data2 = data2[-1].split('GRADIENT_NORM:KCAL/MOL/ANGSTROM')
        data2 = data2[0].replace('=','')
        data2 = data2.replace('D','E')
        data2 = data2.replace('+','')
        data2 = data2.replace('\n','')
        data2 = float(data2)

        HEAT_OF_FORMATION = data2
    except:
        HEAT_OF_FORMATION = None


    #.ATOM_X_OPT
    try:
        data2 = data.split('ATOM_X_OPT:ANGSTROMS')
        data2 = data2[-1].split('ATOM_CHARGES')
        data2 = data2[0] 

        data2 = data2.split('=')
        data2 = data2[-1].split()

        ATOM_X_OPT = data2
    except:
        ATOM_X_OPT = None
        
        
    #,GRADIENTS
    try: 
        data2 = data.split('GRADIENTS:KCAL/MOL/ANGSTROM')
        data2 = data2[-1].split('OVERLAP_MATRIX')
        data2 = data2[0] 
        data2 = data2.split('=')
        data2 = data2[-1].split()

        GRADIENTS = data2
    except:
        GRADIENTS = None
    
    #,charges
    try: 
        data2 = data.split('ATOM_CHARGES')
        data2 = data2[-1].split('AO_CHARGES')
        data2 = data2[0] 
        data2 = data2.split('=')
        data2 = data2[-1].split()
                
        CHARGES = data2
    except:
        CHARGES = None    
    
    #print (HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS)
    return HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS, CHARGES
    
def read_MOPAC_arc (filein = None):
    """ Function doc """
    
def read_ORCA_log (filein = None):
    """ Function doc """
    
def read_XTB_log (filein = None):
    """ Function doc """
    

def read_MOL2 (filein):
    """ Function doc """
    #at  = MolecularProperties.AtomTypes()
    
    with open(filein, 'r') as mol2_file:
        filetext = mol2_file.read()

        molecules     =  filetext.split('@<TRIPOS>MOLECULE')
        firstmolecule =  molecules[1].split('@<TRIPOS>ATOM')
        header        =  firstmolecule[0]
        firstmolecule =  firstmolecule[1].split('@<TRIPOS>BOND')
        raw_atoms     =  firstmolecule[0]
        bonds         =  firstmolecule[1]


    header    = header.split('\n')
    raw_atoms = raw_atoms.split('\n')
    bonds     = bonds.split('\n')
    
    atoms = {
            'id'     : [],
            'name'   : [],
            'xyz'    : [],
            'type'   : [],
            'resi'   : [],
            'resn'   : [],
            'charges': [],
            }
    
    _bonds = []
    
    #pprint (raw_atoms)
    for atom in raw_atoms[1:]:
        if len(atom) == 0:
            pass
        else:
            try:
                a = atom.split()
                atoms['id'     ].append(int(a[0]))
                atoms['name'   ].append(a[1])
                atoms['xyz'    ].append([ float(a[2]), float(a[3]), float(a[4])] )
                atoms['type'   ].append(a[5])
                atoms['resi'   ].append(int(a[6]))
                atoms['resn'   ].append(a[7])
                atoms['charges'].append(float(a[8]))
            except:
                pass
    #pprint (atoms)
    
    
    for bond in bonds:
        b = bond.split()
        
        if len(b) == 4:
            _bonds.append(b)
            
        else:
            pass
    #pprint (_bonds)
    return atoms, _bonds

#read_MOL2 (filein = "/home/fernando/Desktop/NoName.mol2")
