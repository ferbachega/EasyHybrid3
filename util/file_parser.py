from pprint import pprint


def get_file_type(filename):
    file_type = filename.split('.')
    return file_type[-1]
    
    

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
            a = atom.split()
            atoms['id'     ].append(int(a[0]))
            atoms['name'   ].append(a[1])
            atoms['xyz'    ].append([float(a[2]),float(a[3]),float(a[4])])
            atoms['type'   ].append(a[5])
            atoms['resi'   ].append(int(a[6]))
            atoms['resn'   ].append(a[7])
            atoms['charges'].append(float(a[8]))
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
