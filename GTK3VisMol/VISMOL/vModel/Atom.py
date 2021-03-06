#import VISMOL.vModel.atom_types as at 

#GTK3EasyMol/VISMOL/Model/atom_types.py
class Atom:
    """ Class doc """
    
    def __init__ (self, name         ='Xx',
                        index        =None, 
                        symbol       ='X', 
                        pos          = None, 
                        resi         = None, 
                        resn         = None, 
                        chain        = ''  , 
                        atom_id      = 0   , 
                        residue      = None,
                        
                        occupancy    = 0.0,
                        bfactor      = 0.0, 
                        charge       = 0.0,
						bonds_indices= [] ,
                        Vobject      = None):
 
        if pos is None:
            pos = np.array([0.0, 0.0, 0.0])
		
        
        self.pos        = pos     # - coordinates of the first frame
        self.index      = index   # - Remember that the "index" attribute refers to the numbering of atoms 
                                  # (it is not a zero base, it starts at 1 for the first atom)
					    
        self.name       = name    #
        self.symbol     = symbol  #
        self.resi       = resi    #
        self.resn       = resn    #
        self.chain      = chain   #
        self.Vobject    = Vobject
        self.residue    = residue    
        
        
        self.occupancy  = occupancy
        self.bfactor    = bfactor  
        self.charge     = charge   
        
        
        at = Vobject.vismol_session.vConfig.atom_types

        self.atom_id = atom_id        # An unique number
        self.color   = at.get_color    (name)
        self.color_id = None
        #self.color.append(0.0)  

        self.col_rgb      = at.get_color_rgb(name)
        self.radius       = at.get_radius   (name)
        self.vdw_rad      = at.get_vdw_rad  (name)
        self.cov_rad      = at.get_cov_rad  (name)
        self.ball_radius  = at.get_ball_rad (name)

        self.lines          = True
        self.dots           = False
        self.nonbonded      = False
        self.ribbons        = False
        self.ball_and_stick = False
        self.sticks         = False
        self.spheres        = False
        self.surface        = False
        self.bonds_indices  = bonds_indices
        self.bonds          = []
    
    def coords (self, frame = None):
        """ Function doc """
        if frame is None:
            frame  = self.Vobject.vismol_session.glwidget.vm_widget.frame
            #frame  = self.Vobject.vismol_session.glwidget.vm_widget._safe_frame_exchange(self.Vobject)
        
        
        
        #if self.Vobject.vismol_session.glwidget.vm_widget.frame <  0:
        #    self.Vobject.vismol_session.glwidget.vm_widget.frame = 0
        #else:
        #    pass
        #
        #if self.Vobject.vismol_session.glwidget.vm_widget.frame >= (len (self.Vobject.frames)-1):
        #    frame = self.Vobject.frames[len (self.Vobject.frames)-2]
        #else:
        #    frame = self.Vobject.frames[self.Vobject.vismol_session.glwidget.vm_widget.frame]
        
        
        coords = [self.Vobject.frames[frame][(self.index-1)*3  ],
        self.Vobject.frames[frame][(self.index-1)*3+1],
        self.Vobject.frames[frame][(self.index-1)*3+2],]

        return coords

    def define_atom_symbol (self, name):
        """ Function doc """
        self.symbol_list = [ 'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne','Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Xx', 'X' ]
        
        if name in self.symbol_list:
            self.symbol = name
        
        
        else:
            newSymbol = ''
            
            if len(name) > 1:             
            
                if name[1] in  ['0','1','2','3','4','5','6','7','8','9']:
                    
                    if name[0].islower():
                        newSymbol = name[0].upper()
                        if newSymbol in self.symbol_list:
                            self.symbol = newSymbol
                        else:
                            self.symbol = 'Xx'
                    
                    else:
                        newSymbol = name[0].upper()
                        if newSymbol in self.symbol_list:
                            self.symbol = newSymbol
                        else:
                            self.symbol = 'Xx'
                
                
                
                
                else:
                    
                    if name[0].islower():                        
                        newSymbol += name[0].upper()
    
         
            else:
                
                if name[0].islower():
                    newSymbol = name[0].upper()
                    if newSymbol in self.symbol_list:
                        self.symbol = newSymbol
                    else:
                        self.symbol = 'Xx'

 
 
 
