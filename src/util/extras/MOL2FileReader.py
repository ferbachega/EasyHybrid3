#===================================================================================================================================
# . Classes and functions to read MOL2 files.
#===================================================================================================================================

from  pCore                 import logFile                  , \
                                   LogFileActive            , \
                                   TextFileReader
from  pMolecule             import Atom                     , \
                                   Bond                     , \
                                   BondType                 , \
                                   Connectivity             , \
                                   ConvertInputConnectivity , \
                                   Sequence                 , \
                                   System
from  pScientific           import PeriodicTable
from  pScientific.Arrays    import Array
from  pScientific.Geometry3 import Coordinates3
from  pScientific.Symmetry  import PeriodicBoundaryConditions         , \
                                   SpaceGroup_CrystalSystemFromNumber , \
                                   SymmetryParameters
from .ExportImport          import _Importer

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Bond type definitions.
_MOL2BondTypes = { "1"  : ( BondType.Single    , False ) ,
                   "2"  : ( BondType.Double    , False ) ,
                   "3"  : ( BondType.Triple    , False ) ,
                   "am" : ( BondType.Single    , False ) ,
                   "ar" : ( BondType.Undefined , True  ) ,
                   "du" : ( BondType.Undefined , False ) ,
                   "un" : ( BondType.Undefined , False ) ,
                   "nc" : ( None               , False ) }

# . Generic (non-element) atom types.
_GenericAtomTypes = ( "Any", "Du", "Hal", "Het", "Hev", "LP" )

# . Section header.
_RTIString  = "@<TRIPOS>"

#===================================================================================================================================
# . MOL2 file reader class.
#===================================================================================================================================
class MOL2FileReader ( TextFileReader ):
    """MOL2FileReader is the class for MOL2 files that are to be read."""

    _classLabel = "MOL2 File Reader"

    def AtomicNumberFromAtomType ( self, atomType ):
        """Get an atomic number from atom type."""
        atomicNumber = -1
        token        = atomType.split ( ".", 1 )[0]
        if token not in _GenericAtomTypes: atomicNumber = PeriodicTable.AtomicNumber ( token )
        return atomicNumber

    def GetLine ( self, signalWarnings = False ):
        """Get a line."""
        try:
            while True:
                line = next ( self.file ).strip ( )
                self.linesParsed += 1
                if ( len ( line ) > 0 ) and ( not line.startswith ( "#" ) ): break
            return line
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            # . Parse the data.
            try:
                # . Parse all entries.
                while True:
                    line = self.GetLine ( )
                    if   line.startswith ( _RTIString + "ATOM"     ): self.ParseAtomSection     ( )
                    elif line.startswith ( _RTIString + "BOND"     ): self.ParseBondSection     ( )
                    elif line.startswith ( _RTIString + "CRYSIN"   ): self.ParseCrysinSection   ( )
                    elif line.startswith ( _RTIString + "MOLECULE" ): self.ParseMoleculeSection ( )
            except EOFError:
                pass
            # . Complete the connectivity.
            ConvertInputConnectivity ( self.connectivity, {} )
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    def ParseAtomSection ( self ):
        """Parse the ATOM section."""
        if hasattr ( self, "numberOfAtoms" ):
            self.connectivity = Connectivity ( )
            
            atomicCharges = Array.WithExtent        ( self.numberOfAtoms ) ; atomicCharges.Set ( 0.0 )
            xyz           = Coordinates3.WithExtent ( self.numberOfAtoms ) ; xyz.Set           ( 0.0 )
            
            atom_counter    = 1
            residue_counter = 0     #.add by bachega in Sept 10/2024
            previous_numebr = 0     #.add by bachega in Sept 10/2024
            self.atomPaths = []          #.add by bachega in Sept 10/2024
            self.atomTypes = []
            
            for i in range ( self.numberOfAtoms ):
                items = self.GetTokens ( converters = [ int, None, float, float, float, None, None, None, float ] )
                if len ( items ) < 6:
                    self.Warning ( "Invalid ATOM line.", True )
                else:
                    
                    
                    symbol = self.get_symbol ( atom_name = items[1], atom_type = items[5])
                    print(symbol)
                    atomicNumber = PeriodicTable.AtomicNumber ( symbol )
                    #print(atomicNumber, symbol, items[1], items[5])
                    
                    self.atomTypes.append(items[5])
                    #if len(items[5]) > 1:
                    #    atomlabel  = items[5].split('.')
                    #    atomlabel  = atomlabel[0]
                    #else:
                    #    atomlabel  = items[5]
                    
                    #atomicNumber = PeriodicTable.AtomicNumber ( atomlabel )
                    #if atomicNumber <= 0: atomicNumber = self.AtomicNumberFromAtomType ( items[5] )
                    

                    #self.connectivity.AddNode ( Atom.WithOptions ( atomicNumber = atomicNumber, label = items[1] ) )
                    self.connectivity.AddNode ( Atom.WithOptions ( atomicNumber = atomicNumber, label = items[1] ) )
                    #print(atomicNumber)
                    xyz[i,0] = items[2]
                    xyz[i,1] = items[3]
                    xyz[i,2] = items[4]
                    
                    #--------------------------------------------------------
                    #           add by bachega in Sept 10/2024
                    #
                    resi_raw = int(items[6])
                    atomPath = ':'
                    
                    if resi_raw != previous_numebr:
                        residue_counter += 1
                        previous_numebr = resi_raw
                        atom_counter = 1
                        
                    #if len(items[1]) == 1:
                    #    #items[1] =+ str(atom_counter)
                    #    atomlabel = items[1]+str(atom_counter)
                    #else:
                    #    atomlabel = items[1]
                        
                    atomPath += items[7][:3]+'.'+str(residue_counter)+':'+items[1]
                    self.atomPaths.append(atomPath)
                    atom_counter += 1
                    #--------------------------------------------------------
                    #print(items[5], items[6], items[7], atomPath)
                    #print(items[0], atomPath, items[5], items[6])
                    
                    if len ( items ) >= 9: atomicCharges[i] = items[8]
            
            #print(self.atomPaths, len(self.atomPaths))
            #self.sequence = Sequence.FromAtomPaths ( atomPaths ) # add by bachega in Sept 10/2024
            #print(self.sequence)
            
            self.atomicCharges = atomicCharges
            self.xyz           = xyz
            #self.atomTypes    = atom_types
            #print(list(self.atomTypes))
        else:
            self.Warning ( "Unknown number of atoms in molecule.", True )

    def ParseBondSection ( self ):
        """Parse the BOND section."""
        if hasattr ( self, "numberOfBonds" ) and hasattr ( self, "connectivity" ):
            numberOfAtoms = getattr ( self, "numberOfAtoms", -1 )
            for i in range ( self.numberOfBonds ):
                items = self.GetTokens ( converters = [ int, int, int, None ] )
                if len ( items ) < 4:
                    self.Warning ( "Invalid BOND line.", True )
                else:
                    atom1 = items[1] - 1
                    atom2 = items[2] - 1
                    ( bondType, bondIsAromatic ) = _MOL2BondTypes.get ( items[3], ( BondType.Undefined, False ) )
                    #print(bondType, bondIsAromatic)
                    if ( atom1 < 0 ) or ( atom1 >= self.numberOfAtoms ) or ( atom2 < 0 ) or ( atom2 >= self.numberOfAtoms ):
                        self.Warning ( "Bond atom indices out of range: {:d}, {:d}.".format ( atom1, atom2 ), True )
                    if bondType is not None:
                        self.connectivity.AddEdge ( Bond.WithNodes ( self.connectivity.nodes[atom1], self.connectivity.nodes[atom2], isAromatic = bondIsAromatic, type = bondType ) )
        else:
            self.Warning ( "Unknown number of bonds in molecule.", True )

    def ParseCrysinSection ( self ):
        """Parse the CRYSIN section."""
        # . Items are a, b, c, alpha, beta, gamma, space group number and crystal setting.
        items  = self.GetTokens ( converters = 6 * [ float ] + 2 * [ int ] )
        try:
            setting    = items.pop ( -1 )
            spaceGroup = items.pop ( -1 )
            self.crystalSystem      = SpaceGroup_CrystalSystemFromNumber ( spaceGroup )
            self.symmetryParameters = SymmetryParameters ( )
            self.symmetryParameters.SetCrystalParameters ( *items )
        except:
            self.Warning ( "Invalid CRYSIN line.", False )

    def ParseMoleculeSection ( self ):
        """Parse the MOLECULE section."""
        # . Just parse the first two lines for the moment.
        self.label = self.GetLine ( )                             # . Molecule name.
        items      = self.GetTokens ( converters = [ int, int ] ) # . Number of atoms and bonds.
        if len ( items ) > 0: self.numberOfAtoms = items[0]
        if len ( items ) > 1: self.numberOfBonds = items[1]

    @classmethod
    def PathToAtomNames ( selfClass, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToAtomNames ( )

    @classmethod
    def PathToCharges ( selfClass, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToCharges ( )

    @classmethod
    def PathToCoordinates3 ( selfClass, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToCoordinates3 ( )

    @classmethod
    def PathToSystem ( selfClass, path, log = logFile ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToSystem ( )

    def ToAtomNames ( self ):
        """Return atom names."""
        atomNames = None
        if self.isParsed and hasattr ( self, "connectivity" ):
            atomNames = [ atom.label for atom in self.connectivity.nodes ]
        return atomNames

    def ToCharges ( self ):
        """Return charges."""
        charges = None
        if self.isParsed: return getattr ( self, "atomicCharges", None )
        return charges

    def ToCoordinates3 ( self ):
        """Return a coordinates3 object."""
        xyz = None
        if self.isParsed: return getattr ( self, "xyz", None )
        return xyz

    def ToSystem ( self ):
        """Return a system."""
        system = None
        if self.isParsed:
            #try:
            #if hasattr ( self, "atomPaths" ): sequence = Sequence.FromAtomPaths ( atomPaths, atoms = self.atoms )
            
            
            '''
            #.Original from Martin 
            # Needs:
            #  -self.atoms (pDynamo atom objects)
            #  -self.bonds (pDynamo bond objects)
            
            
            system = None
            if self.isParsed:
                
                #if hasattr ( self, "atomPaths" ):  
                #    print(self.atomPaths)
                #else:
                #    print('self.atomPaths is not defined')
                
                try:
                    if hasattr ( self, "atomPaths" ): 
                        sequence = Sequence.FromAtomPaths ( atomPaths, atoms = self.atoms )
                    else:                             
                        sequence = Sequence.FromAtoms ( self.atoms )
                    
                    #print(atomPaths)
                    system = System.FromSequence ( sequence, bonds = ( self.bonds if hasattr ( self, "bonds" ) else None ) )
                    
                    system.label        = self.label
                    system.coordinates3 = self.xyz
                    if hasattr ( self, "crystalSystem" ) and hasattr ( self, "symmetryParameters" ):
                        parameters                = self.crystalSystem.GetUniqueSymmetryParameters ( self.symmetryParameters )
                        system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( self.crystalSystem )
                        system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( **parameters )
                except Exception as e:
                    raise TextFileReaderError ( "Error generating system from MOL2 file.", e.args )
            #'''
            
            
            
            #'''
            try:
                '''
                The function Sequence.FromAtomPaths also requires the list   
                of atoms generated by pDynamo (and subsequently, the bonds). 
                It has not been working very well.
                
                I will keep this parsing system instead of the original one 
                — contact Martin later.
                '''
                #-------------------------------------------------------
                # this is the original one
                #else:
                #    print('self.atomPaths is not defined')
                #try:
                #    if hasattr ( self, "atomPaths" ): sequence = Sequence.FromAtomPaths ( atomPaths, atoms = self.atoms )
                #    else:                             sequence = Sequence.FromAtoms ( self.atoms )
                #    #print(atomPaths)
                #    system = System.FromSequence ( sequence, bonds = ( self.bonds if hasattr ( self, "bonds" ) else None ) )
                #    system.label        = self.label
                #    system.coordinates3 = self.xyz
                
                #-------------------------------------------------------
                
                system              = System.FromConnectivity ( connectivity= self.connectivity)#, withSequence =True)
                system.label        = self.label
                system.coordinates3 = self.xyz
                
                if hasattr ( self, "crystalSystem" ) and hasattr ( self, "symmetryParameters" ):
                    parameters                = self.crystalSystem.GetUniqueSymmetryParameters ( self.symmetryParameters )
                    system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( self.crystalSystem )
                    system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( **parameters )

            except:
                pass
            #'''
        

        '''
        sequence = Sequence.FromAtoms ( system.atoms )
        print(sequence)
        for atom in system.atoms:
            resName, resSeq, iCode = sequence.ParseLabel ( atom.parent.label, fields = 3 )
            print(resName, resSeq, iCode, atom.label)
        '''
        #.Bachega
        # These two attributes are important for EasyHybrid to be able 
        # to discern the sequence in systems generated from connectivity.
        system.sequence_from_mol2 = self.atomPaths
        system.atomType_FromMol2  = self.atomTypes
        
        return system

    def get_symbol (self, atom_name = None, atom_type = None):
        """ Function doc """
        #print(atom_name , atom_type )
        
        #.If there is only one character, this is the chemical symbol itself.
        #. O, H, C, P, S and so...
        if len(atom_name) > 1:
            #.If the first character is numeric, we have an error (anomalous notation)
            if atom_name[0].isalpha:
                #.If the second character is numeric, then the first character is the chemical symbol itself.
                #.like C1, H1, O3  and so... 
                if atom_name[1].isalpha:
                    
                    #.Second atom is upper? like CE2 CA CD1 - you need to check
                    if atom_name[1].isupper:
                        
                        #. In this case we have to check not only the 
                        #  name of the atom, but also the atomic type, 
                        #  because there is a duplicity, for example the 
                        #  name of the atom CD can be confused between 
                        #  cadmium and carbon delta.
                        
                        atom_type = atom_type.split('.') 
                        #. When the topology is written by Avogadro, for example, 
                        #  the atomic type is written according to the hybridization type.
                        #. like O.3 , O.co2 , N.am , C.ar , N.pl3
                        if len (atom_type) == 2:
                            atom_type = atom_type[0] 
                            symbol = atom_type
                        
                        #.When the topology is written by gabedit, for 
                        # example, the atomic type is written using 
                        # AMBER /CHARMM like notation.
                        # CT NA CC C*
                        else:
                            
                            atom_type = atom_type[0] # <-- kind of ['CD'] or ['C','3']
                            #print(atom_type, atom_name)
                            
                            #.If there is only one character, this is 
                            # the chemical symbol itself.
                            # CT, CN, NA, CZ, OH and so on.
                            if len (atom_type) > 1:
                                
                                # CT, CN, NA, CZ, OH and so on.
                                if atom_type[1].isalpha():
                                    
                                    #.If the second character is the same 
                                    # for the atom name and the atomic type, 
                                    # regardless of case, then this is the symbol itself.
                                    if atom_type[1].lower() == atom_name[1].lower():
                                        #.like atom_type = Cd 
                                        #      atom_name = CD
                                        # is cadmium.  
                                        symbol = atom_type[0]+atom_name[1].lower()
                                    
                                    else:
                                        #.like atom_type = CT 
                                        #      atom_name = CD
                                        # is carbom.  
                                        symbol = atom_type[0]
                                
                                
                                else:
                                    #.like when atom_type equal H1, B1, C3, N3, O1, S1, P2...
                                    symbol = atom_type[0]
                            
                            # like when atom_type equal H, B, C, N, O, S, P...
                            else:
                                symbol = atom_type[0]
                            
                    # when atom_name is like Mg, Ca, Cd, Fe, Co and so on... 
                    else:
                        symbol = atom_name
                
                # like H01, 
                else:
                    symbol = atom_name[0]
            
            #.Error atom name connot start with a number 
            else:
                symbol = None    
        
        # like when atom_name equal H, B, C, N, O, S, P...
        else:
            symbol = atom_name
        return symbol   
        


#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : MOL2FileReader.PathToCoordinates3 ,
                         System       : MOL2FileReader.PathToSystem       } ,
                       [ "mol2", "MOL2" ], "Tripos MOL2", defaultFunction = MOL2FileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
