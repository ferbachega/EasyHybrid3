"""The MOPAC QC model."""

import glob, math, os, os.path, subprocess, re

from  pCore                     import logFile           , \
                                       LogFileActive     , \
                                       NotInstalledError
from  pScientific               import PeriodicTable     , \
                                       Units
from  pScientific.Arrays        import Array
from  pScientific.Geometry3     import Coordinates3      , \
                                       Vector3
from  pScientific.RandomNumbers import RandomString
from .QCDefinitions             import ChargeModel
from .QCModel                   import QCModel           , \
                                       QCModelState
from .QCModelError              import QCModelError

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Default error suffix.
_DefaultErrorPrefix = "error_"

# . Default job name.
_DefaultJobName = "MOPACJob"

# . Command environment variable.
_MOPACCommand = "PDYNAMO3_MOPACCOMMAND"

# . Scratch directory.
_MOPACScratch = os.path.join ( os.getenv ( "PDYNAMO3_SCRATCH" ), "MOPACScratch" )

HARTREE = 627.509474 #kcal mol-1
BOHR    = 0.529177   # Angtrons

KCAL_TO_HARTREE = 1/HARTREE # from kcal to hartree
GRAD_FACTOR     = (KCAL_TO_HARTREE) / (1/BOHR) # GRAD factor to convert from KCAL/ANGSTROM  to Hartree/Borh



#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelMOPACState ( QCModelState ):
    """A QC model state."""

    _attributable = dict ( QCModelState._attributable )
    _attributable.update ( { "deleteJobFiles" : True , 
                             "paths"          : None } )

    def __del__ ( self ):
        """Deallocation."""
        self.DeleteJobFiles ( )

    def DeleteJobFiles ( self ):
        """Delete job files."""
        if self.deleteJobFiles:
            try:
                jobFiles = glob.glob ( os.path.join ( self.paths["Glob"] + "*" ) )
                for jobFile in jobFiles: os.remove ( jobFile )
                scratch  = self.paths.get ( "Scratch", None )
                if scratch is not None: os.rmdir ( scratch ) # . Only deleted if random.
            except:
                pass

    def DeterminePaths ( self, scratch, deleteJobFiles = True, randomJob = False, randomScratch = False ):
        """Determine the paths needed by an MOPAC job."""
        paths = {}
        if randomJob: job = RandomString ( )
        else:         job = _DefaultJobName
        if randomScratch:
            scratch          = os.path.join ( scratch, RandomString ( ) )
            paths["Scratch"] = scratch # . Only set if random.
        if not os.path.exists ( scratch ): os.mkdir ( scratch )
        jobRoot       = os.path.join ( scratch, job )
        paths["Glob"] = jobRoot
        for ( key, ext ) in ( #( "EnGrad" , "engrad" ) ,
                              ( "Input"  , "mop"    ) ,
                              ( "Output" , "aux"    ) ,
                              ( "Output2" , "out"    ) ,
                              #( "PCGrad" , "pcgrad" ) ,
                              #( "PC"     , "pc"     ) ,
                              #( "Coord"  , "coord"  ) 
                              ):
            paths[key] = "{:s}.{:s}".format ( jobRoot, ext )
        # . Finish up.
        self.deleteJobFiles = deleteJobFiles
        self.paths          = paths

    def SaveErrorFiles ( self, message ):
        """Save the input and output files for inspection if there is an error."""
        for key in (  "Input" , "Output"   ):
            path = self.paths[key]
            if os.path.exists ( path ):
                ( head, tail ) = os.path.split ( path )
                os.rename ( path, os.path.join ( head, _DefaultErrorPrefix + tail ) )
        ( head, tail ) = os.path.split ( self.paths["Glob"] )
        raise QCModelError ( message + "\nCheck the files \"{:s}*\".".format ( os.path.join ( head, _DefaultErrorPrefix + tail ) ) )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelMOPAC ( QCModel ):
    """The MOPAC QC model class."""

    _attributable = dict ( QCModel._attributable )
    _classLabel   = "MOPAC QC Model"
    _stateObject  = QCModelMOPACState
    _summarizable = dict ( QCModel._summarizable )
    
    _attributable.update ( { "deleteJobFiles" : False          ,
                             "keywords"       : None           ,
                                                               
                             "method"         : 'PM7'          ,
                             #"parallel"       : 1              ,
                             #"fermi_temp"     : 300            ,
                             #"iterations"     : 300            , 
                             #"vfukui"         : False          ,
                             #"acc"            : 1.0            ,
                                                               
                             "randomJob"      : False          ,
                             "randomScratch"  : False          ,
                             #"lmo"            : False          ,
                             #"json"           : False          ,

                             "scratch"        : _MOPACScratch    } )
    
    
    _summarizable.update ( { "deleteJobFiles" : "Delete Job Files"        ,
                             "method"         : "MOPAC method"           ,
                             #"parallel"       : "Number of CPUs"          ,
                             #"fermi_temp"     : "Fermi Temperature"       , 
                             #"iterations"     : "Number of Iterations"    ,
                             #"vfukui"         : 'Fukui Indices'           ,
                             "scratch"        : 'scratch'                 ,
                             #"json"           : 'Write JSON Logfile'      ,              
                             #"acc"            : 'Accuracy for SCC'        ,
                             #"lmo"            : 'Localization of Orbitals',                             
                             
                             
                             
                             "randomJob"      : "Random Job"         ,
                             "randomScratch"  : "Random Scratch"     } )
        
    def AtomicCharges ( self, target, chargeModel = ChargeModel.Mulliken ):
        """Atomic charges."""
        source = target.scratch.MOPACOutputData
        if   chargeModel is ChargeModel.CHelpG : return source.get ( "Chelpg Charges"   , None )
        elif chargeModel is ChargeModel.Loewdin: return source.get ( "Loewdin Charges"  , None )
        else:                                    return source.get ( "Mulliken Charges" , None )

    def AtomicSpins ( self, target, chargeModel = "Mulliken" ):
        """Atomic spins."""
        source = target.scratch.MOPACOutputData
        if chargeModel is ChargeModel.Loewdin: return source.get ( "Loewdin Spins"  , None )
        else:                                  return source.get ( "Mulliken Spins" , None )

    def BondOrders ( self, target, chargeModel = None ):
        """Bond Orders - Mayer only."""
        return target.scratch.MOPACOutputData.get ( "Mayer Bond Orders", None )

    def BuildModel ( self, target, qcSelection = None ):
        """Build the model."""
        state = super ( QCModelMOPAC, self ).BuildModel ( target, qcSelection = qcSelection )
        state.DeterminePaths ( self.scratch                         ,
                               deleteJobFiles = self.deleteJobFiles ,
                               randomJob      = self.randomJob      ,
                               randomScratch  = self.randomScratch  )
        return state

    def DipoleMoment ( self, target, center = None ):
        """Dipole Moment."""
        return target.scratch.MOPACOutputData.get ( "Dipole", None )

    def Energy ( self, target ):
        """Calculate the quantum chemical energy."""
        doGradients    = target.scratch.doGradients
        MOPACOutputData = {}
        state          = getattr ( target, self.__class__._stateName )
        target.scratch.MOPACOutputData = MOPACOutputData
        
        self.WriteInputFile ( target, doGradients, ( target.nbModel is not None ), target.scratch.qcCoordinates3AU )
        
        isOK = self.Execute ( state, target )
        if not isOK: state.SaveErrorFiles ( "Error executing program." )
        
        if doGradients:
            isOK = self.ReadEngradFile ( target, MOPACOutputData, target.scratch.qcGradients3AU )
            
            if not isOK: state.SaveErrorFiles ( "Error reading engrad file." )
        
        isOK = self.ReadOutputFile ( target, MOPACOutputData ) # . Returns whether converged or an error.
        if not isOK: state.SaveErrorFiles ( "Error reading output file." )
        target.scratch.energyTerms["MOPAC QC"] = ( MOPACOutputData["Energy"] * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )

    def multiplicity_to_unpaired (self, multiplicity = 1):
        """ Function doc """
        unpaired = multiplicity-1
        return unpaired
    
    def Execute ( self, state, target ):
        """Execute the MOPAC job."""
        #print(self._attributable['cpus'])
        #try:
            #outFile = open ( state.paths["Output"], "w" )
            
        backup = os.getcwd()
        directory = os.path.dirname(state.paths["Input"])
        os.chdir(directory)
        outFile = state.paths["Output"]
        #print("Output", outFile)
        # starting MOPAC exec
        cmd = self.command
        
        
        # charge and multiplicity
        charge       = target.electronicState.charge
        multiplicity = target.electronicState.multiplicity
        unpaired = self.multiplicity_to_unpaired(multiplicity)


        #.gradients
        cmd +=  ' '+state.paths["Input"] + ' > /dev/null 2>&1'
        #print('cmd', cmd)
        
        os.system(cmd)
        return True
                
    def OrbitalEnergies ( self, target ):
        """Orbital energies and HOMO and LUMO indices."""
        return ( target.scratch.MOPACOutputData.get ( "Orbital Energies", None ) ,
                 target.scratch.MOPACOutputData.get ( "HOMO"            , -1   ) ,
                 target.scratch.MOPACOutputData.get ( "LUMO"            , -1   ) )

    def ReadEngradFile ( self, target, MOPACOutputData, gradients3 ):
        """Read an engrad file.
        KCAL_TO_HARTREE 
        GRAD_FACTOR     
        """
        # . The energy and gradients are in atomic units.
        state = getattr ( target, self.__class__._stateName )
        #HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS = self.ReadAuxFile (state.paths["Output"])
        HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS = self.ReadOutFile (state.paths["Output2"])
        
        
        #1 kcal/mol to Hartree (Eh): 1 kcal/mol=0.00159362 Eh1 kcal/mol=0.00159362 Eh
        MOPACOutputData["Energy"] = HEAT_OF_FORMATION*KCAL_TO_HARTREE 
        #print(MOPACOutputData["Energy"])
        #                                             0.000843297
        for i in range ( len ( state.atomicNumbers ) ): 
            gradients3[i,0] = float(GRADIENTS[i*3+0])*GRAD_FACTOR #Force in kcal/A˚ to Eh/bohr
            gradients3[i,1] = float(GRADIENTS[i*3+1])*GRAD_FACTOR #Force in kcal/A˚ to Eh/bohr
            gradients3[i,2] = float(GRADIENTS[i*3+2])*GRAD_FACTOR #Force in kcal/A˚ to Eh/bohr
            #gradients3[i,0] = float(GRADIENTS[i*3+0])*0.00084372 #Force in kcal/A˚ to Eh/bohr
            #gradients3[i,1] = float(GRADIENTS[i*3+1])*0.00084372 #Force in kcal/A˚ to Eh/bohr
            #gradients3[i,2] = float(GRADIENTS[i*3+2])*0.00084372 #Force in kcal/A˚ to Eh/bohr
        return True
        
    def ReadOutputFile (self, target, MOPACOutputData):
        """ Function doc 
        KCAL_TO_HARTREE
        """
        state  = getattr ( target, self.__class__._stateName )
        #print(state.paths["Output"])
        scratch         = { "Is Successful" : False }
        #HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS = self.ReadAuxFile ( state.paths["Output"])
        
        HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS = self.ReadOutFile ( state.paths["Output2"])
        
        scratch["Energy"]       = HEAT_OF_FORMATION*KCAL_TO_HARTREE
        scratch         = { "Is Successful" : False }
        MOPACOutputData["Energy"] = HEAT_OF_FORMATION*KCAL_TO_HARTREE
        return True

    def ReadOutFile (self, _file):
        """ Function doc """
        data = open(_file, 'r')
        data = data.readlines()
        
        HEAT_OF_FORMATION = None 
        ATOM_CORE = None
        ATOM_X_OPT = None
        GRADIENTS = []
        
        for line in data:
            # 0     1   2    3       4           5
            #FINAL HEAT OF FORMATION =       -849.59445 KCAL/MOL =   -3554.70318 KJ/MOL
            if 'FINAL HEAT OF FORMATION' in line:
                line2 = line.split()
                HEAT_OF_FORMATION = float(line2[5])
            
            if 'KCAL/ANGSTROM' in line:
                line2 = line.split()
                grad = line2[-2]
                GRADIENTS.append(grad)
        #print(HEAT_OF_FORMATION, GRADIENTS )
        return HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS
    
    def ReadAuxFile (self, auxfile):
        """ Function doc """
        data = open(auxfile, 'r')
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
        
        #print (HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS)
        return HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS
       
    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCModelMOPAC, self ).SummaryItems ( )
        if self.keywords is not None:
            n = len ( self.keywords )
            items.append ( ( "Keywords", "{:s}".format ( "/".join ( self.keywords[0:min(2,n)] ) ) ) )
        return items

    # . Coordinates are written in atomic units (keywords must contain "Bohrs").
    def WriteInputFile ( self, target, doGradients, doQCMM, coordinates3 ):
        """Write an input file."""
        state  = getattr ( target, self.__class__._stateName )
        #print(list(coordinates3[0]))
        '''
        * ===============================
        * Input file for Mopac
        * ===============================
        PM7 1SCF GRADIENTS EPS=78.39 RSOLV=1.3 CHARGE=0 Singlet  BONDS AUX 

        Mopac file generated by Gabedit
        '''
        charge       = target.electronicState.charge
        multiplicity = target.electronicState.multiplicity
        #print(self.method, 'self.method')
        multip = {
                  1:'Singlet',
                  2:'Doublet', 
                  3:'Triblet', 
                  4:'QUARTET',
                  5:'QUINTET'
                  }
        
        
        text = ''
        text +=  '* ===============================\n'
        text +=  '*     MOPAC input file           \n'
        text +=  '* ===============================\n'
        #text +=  '{} 1SCF GRADIENTS CHARGE={} {}  BONDS AUX DENOUT OLDENS\n'.format(self.method, charge, multip[multiplicity])
        text +=  '{} CHARGE={} {} \n'.format(self.method, charge, multip[multiplicity])
        text +=  '\n'
        text +=  'Generated by EasyHybrid\n'
        
        
        
        
        coordFile = open ( state.paths["Input"], "w" )
        coordFile.write ( text)
        for ( i, n ) in enumerate ( state.atomicNumbers ):
            
            #.pDynamo exports the coordinates3 obj in  bohr units 
            #. 1.8897259886 is the factor to convert 1 A to 1 Bohr
            
            #                  .N        X           Y           Z   
            coordFile.write ( "{:<3s} {:20.10f} 1 {:20.10f} 1 {:20.10f} 1\n".format (PeriodicTable.Symbol ( n ) ,
                                                                               coordinates3[i,0]/(1/BOHR) ,
                                                                               coordinates3[i,1]/(1/BOHR) ,
                                                                               coordinates3[i,2]/(1/BOHR) ,
                                                                               #coordinates3[i,0]/1.8897259886 ,
                                                                               #coordinates3[i,1]/1.8897259886 ,
                                                                               #coordinates3[i,2]/1.8897259886 ,
                                                                                )
                                                                            )
                                                                            
            #coordFile.write ( "{:20.10f}{:20.10f}{:20.10f} {:<12s}\n".format ( 
            #                                                                coordinates3[i,0] ,
            #                                                                coordinates3[i,1] ,
            #                                                                coordinates3[i,2] ,
            #                                                                PeriodicTable.Symbol ( n ))
            #                                                                )

        #coordFile.write ( "$end\n")
        coordFile.close ( )
        
    @property
    def command ( self ):
        """Get the command to execute the program."""
        command = self.__dict__.get ( "_command", None )
        if command is None:
            command = os.getenv ( _MOPACCommand )
            # . Command must point to an executable file.
            if  ( command is None ) or not ( os.path.isfile ( command ) and os.access ( command, os.X_OK ) ):
                raise NotInstalledError ( "MOPAC executable not found." )
            else:
                self.__dict__["_command"] = command
        return command

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

