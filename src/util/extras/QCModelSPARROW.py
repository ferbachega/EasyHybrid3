"""The SPARROW QC model."""

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
#from pBabel                     import *
#from pBabel                     import *

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Default error suffix.
_DefaultErrorPrefix = "error_"

# . Default job name.
_DefaultJobName = "SPARROWJob"

# . Command environment variable.
_SPARROWCommand = "PDYNAMO3_SPARROWCOMMAND"

# . Scratch directory.
_SPARROWScratch = os.path.join ( os.getenv ( "PDYNAMO3_SCRATCH" ), "SPARROWScratch" )

HARTREE = 627.509474 #kcal mol-1
BOHR    = 0.529177   # Angtrons

KCAL_TO_HARTREE = 1/HARTREE # from kcal to hartree
GRAD_FACTOR     = (KCAL_TO_HARTREE) / (1/BOHR) # GRAD factor to convert from KCAL/ANGSTROM  to Hartree/Borh



#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelSPARROWState ( QCModelState ):
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
        """Determine the paths needed by an SPARROW job."""
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
                              ( "Input"  , "xyz"    ) ,
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
class QCModelSPARROW ( QCModel ):
    """The SPARROW QC model class."""

    _attributable = dict ( QCModel._attributable )
    _classLabel   = "SPARROW QC Model"
    _stateObject  = QCModelSPARROWState
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

                             "scratch"        : _SPARROWScratch    } )
    
    
    _summarizable.update ( { "deleteJobFiles" : "Delete Job Files"        ,
                             "method"         : "PM6"           ,
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
        source = target.scratch.SPARROWOutputData
        if   chargeModel is ChargeModel.CHelpG : return source.get ( "Chelpg Charges"   , None )
        elif chargeModel is ChargeModel.Loewdin: return source.get ( "Loewdin Charges"  , None )
        else:                                    return source.get ( "Mulliken Charges" , None )

    def AtomicSpins ( self, target, chargeModel = "Mulliken" ):
        """Atomic spins."""
        source = target.scratch.SPARROWOutputData
        if chargeModel is ChargeModel.Loewdin: return source.get ( "Loewdin Spins"  , None )
        else:                                  return source.get ( "Mulliken Spins" , None )

    def BondOrders ( self, target, chargeModel = None ):
        """Bond Orders - Mayer only."""
        return target.scratch.SPARROWOutputData.get ( "Mayer Bond Orders", None )

    def BuildModel ( self, target, qcSelection = None ):
        """Build the model."""
        state = super ( QCModelSPARROW, self ).BuildModel ( target, qcSelection = qcSelection )
        state.DeterminePaths ( self.scratch                         ,
                               deleteJobFiles = self.deleteJobFiles ,
                               randomJob      = self.randomJob      ,
                               randomScratch  = self.randomScratch  )
        return state

    def DipoleMoment ( self, target, center = None ):
        """Dipole Moment."""
        return target.scratch.SPARROWOutputData.get ( "Dipole", None )

    def Energy ( self, target ):
        """Calculate the quantum chemical energy."""
        doGradients    = target.scratch.doGradients
        SPARROWOutputData = {}
        state          = getattr ( target, self.__class__._stateName )
        target.scratch.SPARROWOutputData = SPARROWOutputData
        
        self.WriteInputFile ( target, doGradients, ( target.nbModel is not None ), target.scratch.qcCoordinates3AU )
        
        isOK = self.Execute ( state, target )
        if not isOK: state.SaveErrorFiles ( "Error executing program." )
        
        if doGradients:
            isOK = self.ReadEngradFile ( target, SPARROWOutputData, target.scratch.qcGradients3AU )
            
            if not isOK: state.SaveErrorFiles ( "Error reading engrad file." )
        
        isOK = self.ReadOutputFile ( target, SPARROWOutputData ) # . Returns whether converged or an error.
        if not isOK: state.SaveErrorFiles ( "Error reading output file." )
        target.scratch.energyTerms["SPARROW QC"] = ( SPARROWOutputData["Energy"] * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )

    def multiplicity_to_unpaired (self, multiplicity = 1):
        """ Function doc """
        unpaired = multiplicity-1
        return unpaired
    
    def Execute ( self, state, target ):
        """Execute the SPARROW job."""
        """Execute the xtb job."""
        #print(self._attributable['cpus'])
        #try:
            #outFile = open ( state.paths["Output"], "w" )
        #os.system('/home/fernando/programs/sparrow/install/bin/sparrow')
        backup = os.getcwd()
        #directory = os.path.dirname(state.paths["Coord"])
        #os.chdir(directory)
        #outFile = state.paths["Output"]
        #print(self.command)
        # starting xtb exec
        cmd = self.command
        
        
        # charge and multiplicity
        charge       = target.electronicState.charge
        multiplicity = target.electronicState.multiplicity
        c_and_u  = ' -c {} -s {} '.format(charge, multiplicity)            
        cmd += c_and_u
        
        
        #.CPUs
        #cpus = ' -P {} '.format(self.parallel)
        #cmd += cpus
        
        cmd += ' -M {}'.format('PM6')
        cmd +=  ' -G' 
        cmd +=  ' -o'
        #aditional keys?
        if  self.keywords:
            #for key in self.keywords:
            cmd += ' ' + self.keywords + ' '
        
        # adding inputfile
        cmd +=  ' -x ' + state.paths["Input"] + ' > SPARROWJob.out'
        #print(cmd)
        os.system(cmd)
        return True
                
    def OrbitalEnergies ( self, target ):
        """Orbital energies and HOMO and LUMO indices."""
        return ( target.scratch.SPARROWOutputData.get ( "Orbital Energies", None ) ,
                 target.scratch.SPARROWOutputData.get ( "HOMO"            , -1   ) ,
                 target.scratch.SPARROWOutputData.get ( "LUMO"            , -1   ) )


    def ReadEngradFile(self, target, SPARROWOutputData, gradients3):
        """Read gradients.dat file."""

        state = getattr(target, self.__class__._stateName)

        try:
            gradients = []

            with open("gradients.dat", "r") as data:
                for line in data:

                    if "#" in line:
                        continue

                    values = line.split()

                    if len(values) == 3:
                        gradients.extend(float(v) for v in values)
            
            #print(gradients, len(gradients))
            for i in range(len(state.atomicNumbers)):
                for j in range(3):
                    #print(gradients3[i, j])
                    gradients3[i, j] = gradients[i * 3 + j]
            return True

        except Exception:
            return False

    def ReadEngradFile_old ( self, target, SPARROWOutputData, gradients3 ):
        """Read an engrad file.
        KCAL_TO_HARTREE 
        GRAD_FACTOR     
        """
        # . The energy and gradients are in atomic units.
        state = getattr ( target, self.__class__._stateName )



        data = open('gradients.dat', 'r')
        data = data.readlines()
        GRADIENTS = []
        for line in data:
            if '#' in line:
                pass
            else:
                #print(line)
                line2 = line.split()
                if len(line2) ==3:
                    for value in line2:
                        GRADIENTS.append(float(value))
        #print(GRADIENTS)
        print(len ( state.atomicNumbers ))
        for i in range ( len ( state.atomicNumbers ) ): 
            #gradients3[i,0] = float(GRADIENTS[i*3+0])*GRAD_FACTOR #Force in kcal/A˚ to Eh/bohr
            #gradients3[i,1] = float(GRADIENTS[i*3+1])*GRAD_FACTOR #Force in kcal/A˚ to Eh/bohr
            #gradients3[i,2] = float(GRADIENTS[i*3+2])*GRAD_FACTOR #Force in kcal/A˚ to Eh/bohr
            gradients3[i,0] = float(GRADIENTS[i*3+0]) # Eh/bohr
            gradients3[i,1] = float(GRADIENTS[i*3+1]) # Eh/bohr
            gradients3[i,2] = float(GRADIENTS[i*3+2]) # Eh/bohr
        return True
        
    def ReadOutputFile (self, target, SPARROWOutputData):
        """ Function doc 
        """
        state  = getattr ( target, self.__class__._stateName )
        scratch         = { "Is Successful" : False }

        
        data  = open('energy.dat', 'r')
        for line in data:
            if '#' in line:
                pass
            else:
                line2 = line.split()
                if len(line2) ==1:
                    energy = float(line2[0])
        scratch["Energy"]       = energy
        scratch         = { "Is Successful" : False }
        SPARROWOutputData["Energy"] = energy
        return True

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCModelSPARROW, self ).SummaryItems ( )
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
        * Input file for SPARROW
        * ===============================
        SPARROW file generated by Gabedit
        '''
        #ExportSystem ( 'test2.xyz', target)
        BOHR_TO_ANGSTROM = 0.52917721092
        text = ''
        for ( i, n ) in enumerate ( state.atomicNumbers ):
            
            #.pDynamo exports the coordinates3 obj in  bohr units 
            #. 1.8897259886 is the factor to convert 1 A to 1 Bohr
            
            #                  .N        X           Y           Z   
            text += "{:<3s} {:20.10f} {:20.10f} {:20.10f} \n".format (PeriodicTable.Symbol ( n ) ,
                                                              coordinates3[i,0]*BOHR_TO_ANGSTROM ,
                                                              coordinates3[i,1]*BOHR_TO_ANGSTROM ,
                                                              coordinates3[i,2]*BOHR_TO_ANGSTROM ,
                                                               )
        #print(text)
        text1 =   ''
        text1 +=  '{}\n'.format(i+1)
        text1 +=  '\n'
        text  = text1+text
        
        coordFile = open ( state.paths["Input"], "w" )        
        coordFile.write ( text)
        coordFile.close ( )


    #def ReadOutFile (self, _file):
    #    """ Function doc """
    #    data = open(_file, 'r')
    #    data = data.readlines()
    #    
    #    HEAT_OF_FORMATION = None 
    #    ATOM_CORE = None
    #    ATOM_X_OPT = None
    #    GRADIENTS = []
    #    
    #    for line in data:
    #        # 0     1   2    3       4           5
    #        #FINAL HEAT OF FORMATION =       -849.59445 KCAL/MOL =   -3554.70318 KJ/MOL
    #        if 'FINAL HEAT OF FORMATION' in line:
    #            line2 = line.split()
    #            HEAT_OF_FORMATION = float(line2[5])
    #        
    #        if 'KCAL/ANGSTROM' in line:
    #            line2 = line.split()
    #            grad = line2[-2]
    #            GRADIENTS.append(grad)
    #    #print(HEAT_OF_FORMATION, GRADIENTS )
    #    return HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS
    #
    #def ReadAuxFile (self, auxfile):
    #    """ Function doc """
    #    data = open(auxfile, 'r')
    #    data = data.read()
    #    
    #    #.ATOM_CORE
    #    try:
    #        data2 = data.split('ATOM_CORE')
    #        data2 = data2[-1].split('ATOM_X:ANGSTROMS')
    #        data2 = data2[0] 
    #
    #        data2 = data2.split('=')
    #        data2 = data2[-1].split()
    #
    #        ATOM_CORE = data2
    #    except:
    #        ATOM_CORE = None
    #
    #    #.HEAT_OF_FORMATION
    #    try:
    #        data2 = data.split('HEAT_OF_FORMATION:KCAL/MOL')
    #        data2 = data2[-1].split('GRADIENT_NORM:KCAL/MOL/ANGSTROM')
    #        data2 = data2[0].replace('=','')
    #        data2 = data2.replace('D','E')
    #        data2 = data2.replace('+','')
    #        data2 = data2.replace('\n','')
    #        data2 = float(data2)
    #
    #        HEAT_OF_FORMATION = data2
    #    except:
    #        HEAT_OF_FORMATION = None
    #
    #
    #    #.ATOM_X_OPT
    #    try:
    #        data2 = data.split('ATOM_X_OPT:ANGSTROMS')
    #        data2 = data2[-1].split('ATOM_CHARGES')
    #        data2 = data2[0] 
    #
    #        data2 = data2.split('=')
    #        data2 = data2[-1].split()
    #
    #        ATOM_X_OPT = data2
    #    except:
    #        ATOM_X_OPT = None
    #        
    #        
    #    #,GRADIENTS
    #    try: 
    #        data2 = data.split('GRADIENTS:KCAL/MOL/ANGSTROM')
    #        data2 = data2[-1].split('OVERLAP_MATRIX')
    #        data2 = data2[0] 
    #        data2 = data2.split('=')
    #        data2 = data2[-1].split()
    #
    #        GRADIENTS = data2
    #    except:
    #        GRADIENTS = None
    #    
    #    #print (HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS)
    #    return HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS
       

        
    @property
    def command ( self ):
        """Get the command to execute the program."""
        command = self.__dict__.get ( "_command", None )
        if command is None:
            command = os.getenv ( _SPARROWCommand )
            # . Command must point to an executable file.
            if  ( command is None ) or not ( os.path.isfile ( command ) and os.access ( command, os.X_OK ) ):
                raise NotInstalledError ( "SPARROW executable not found." )
            else:
                self.__dict__["_command"] = command
        return command

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

