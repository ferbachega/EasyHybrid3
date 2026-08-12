"""The XTB QC model."""

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
_DefaultJobName = "XTBJob"

# . Command environment variable.
_XTBCommand = "PDYNAMO3_XTBCOMMAND"

# . Scratch directory.
_XTBScratch = os.path.join ( os.getenv ( "PDYNAMO3_SCRATCH" ), "XTBScratch" )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelXTBState ( QCModelState ):
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
        """Determine the paths needed by an XTB job."""
        paths = {}
        if randomJob: job = RandomString ( )
        else:         job = _DefaultJobName
        if randomScratch:
            scratch          = os.path.join ( scratch, RandomString ( ) )
            paths["Scratch"] = scratch # . Only set if random.
        if not os.path.exists ( scratch ): os.mkdir ( scratch )
        jobRoot       = os.path.join ( scratch, job )
        paths["Glob"] = jobRoot
        for ( key, ext ) in ( ( "EnGrad" , "engrad" ) ,
                              ( "Input"  , "inp"    ) ,
                              ( "Output" , "log"    ) ,
                              ( "PCGrad" , "pcgrad" ) ,
                              ( "PC"     , "pc"     ) ,
                              ( "Coord"  , "coord"  ) ):
            paths[key] = "{:s}.{:s}".format ( jobRoot, ext )
        # . Finish up.
        self.deleteJobFiles = deleteJobFiles
        self.paths          = paths

    def SaveErrorFiles ( self, message ):
        """Save the input and output files for inspection if there is an error."""
        for key in ( "EnGrad" , "Input" , "Output" , "PCGrad" , "PC", 'Coord' ):
            path = self.paths[key]
            if os.path.exists ( path ):
                ( head, tail ) = os.path.split ( path )
                os.rename ( path, os.path.join ( head, _DefaultErrorPrefix + tail ) )
        ( head, tail ) = os.path.split ( self.paths["Glob"] )
        raise QCModelError ( message + "\nCheck the files \"{:s}*\".".format ( os.path.join ( head, _DefaultErrorPrefix + tail ) ) )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelXTB ( QCModel ):
    """The XTB QC model class."""

    _attributable = dict ( QCModel._attributable )
    _classLabel   = "XTB QC Model"
    _stateObject  = QCModelXTBState
    _summarizable = dict ( QCModel._summarizable )
    
    _attributable.update ( { "deleteJobFiles" : False          ,
                             "keywords"       : None           ,
                                                               
                             "gfn"            : 2              ,
                             "parallel"       : 1              ,
                             "fermi_temp"     : 300            ,
                             "iterations"     : 300            , 
                             "vfukui"         : False          ,
                             "acc"            : 1.0            ,
                                                               
                             "randomJob"      : False          ,
                             "randomScratch"  : False          ,
                             #"lmo"            : False          ,
                             #"json"           : False          ,

                             "scratch"        : _XTBScratch    } )
    
    
    _summarizable.update ( { "deleteJobFiles" : "Delete Job Files"        ,
                             "gfn"            : "GFNn-xTB Type"           ,
                             "parallel"       : "Number of CPUs"          ,
                             "fermi_temp"     : "Fermi Temperature"       , 
                             "iterations"     : "Number of Iterations"    ,
                             "vfukui"         : 'Fukui Indices'           ,
                             "scratch"        : 'scratch'                 ,
                             #"json"           : 'Write JSON Logfile'      ,              
                             "acc"            : 'Accuracy for SCC'        ,
                             #"lmo"            : 'Localization of Orbitals',                             
                             
                             
                             
                             "randomJob"      : "Random Job"         ,
                             "randomScratch"  : "Random Scratch"     } )
    
    #self.gfn  = _attributable['gfn']
    #self.cpus = _attributable['cpus']

    
    def AtomicCharges ( self, target, chargeModel = ChargeModel.Mulliken ):
        """Atomic charges."""
        source = target.scratch.XTBOutputData
        if   chargeModel is ChargeModel.CHelpG : return source.get ( "Chelpg Charges"   , None )
        elif chargeModel is ChargeModel.Loewdin: return source.get ( "Loewdin Charges"  , None )
        else:                                    return source.get ( "Mulliken Charges" , None )

    def AtomicSpins ( self, target, chargeModel = "Mulliken" ):
        """Atomic spins."""
        source = target.scratch.XTBOutputData
        if chargeModel is ChargeModel.Loewdin: return source.get ( "Loewdin Spins"  , None )
        else:                                  return source.get ( "Mulliken Spins" , None )

    def BondOrders ( self, target, chargeModel = None ):
        """Bond Orders - Mayer only."""
        return target.scratch.XTBOutputData.get ( "Mayer Bond Orders", None )

    def BuildModel ( self, target, qcSelection = None ):
        """Build the model."""
        state = super ( QCModelXTB, self ).BuildModel ( target, qcSelection = qcSelection )
        state.DeterminePaths ( self.scratch                         ,
                               deleteJobFiles = self.deleteJobFiles ,
                               randomJob      = self.randomJob      ,
                               randomScratch  = self.randomScratch  )
        return state

    def DipoleMoment ( self, target, center = None ):
        """Dipole Moment."""
        return target.scratch.XTBOutputData.get ( "Dipole", None )

    def Energy ( self, target ):
        """Calculate the quantum chemical energy."""
        doGradients    = target.scratch.doGradients
        XTBOutputData = {}
        state          = getattr ( target, self.__class__._stateName )
        target.scratch.XTBOutputData = XTBOutputData
        
        self.WriteInputFile ( target, doGradients, ( target.nbModel is not None ), target.scratch.qcCoordinates3AU )
        
        isOK = self.Execute ( state, target )
        if not isOK: state.SaveErrorFiles ( "Error executing program." )
        if doGradients:
            isOK = self.ReadEngradFile ( target, XTBOutputData, target.scratch.qcGradients3AU )
            #print ('\n\n')
            #print (target.scratch.qcGradients3AU)
            #print ('\n\n')
            if not isOK: state.SaveErrorFiles ( "Error reading engrad file." )
        isOK = self.ReadOutputFile ( target, XTBOutputData ) # . Returns whether converged or an error.
        if not isOK: state.SaveErrorFiles ( "Error reading output file." )
        target.scratch.energyTerms["XTB QC"] = ( XTBOutputData["Energy"] * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )

    def multiplicity_to_unpaired (self, multiplicity = 1):
        """ Function doc """
        unpaired = multiplicity-1
        return unpaired
    
    def Execute ( self, state, target ):
        """Execute the xtb job."""
        

        backup = os.getcwd()
        directory = os.path.dirname(state.paths["Coord"])
        os.chdir(directory)
        outFile = state.paths["Output"]



        # O erro é anteriror
        ##Check xTB tmp path.
        #if os.path.isdir(directory):
        #    print("A pasta existe")
        #    #return False
        #else:
        #    print("A pasta não existe")
        #    PDYNAMO3_SCRATCH = os.getenv ( "PDYNAMO3_SCRATCH" )
        #    tag = 'XTB_'+ RandomString ( )
        #    scratch    = os.path.join ( scratch, tag )
        #    
        #    print("Nova pasta:", scratch)
        #    #target.Summary()
        #    target.qcState.DeterminePaths (scratch)
            

        # Build the command as an ARGUMENT LIST (not a shell string). Safer than
        # os.system (no shell parsing) and lets us run with cwd=directory instead
        # of chdir'ing the whole process -- removing the old cwd-leak bug.
        args = [ self.command ]
        
        # charge and multiplicity
        charge       = target.electronicState.charge
        multiplicity = target.electronicState.multiplicity
        unpaired = self.multiplicity_to_unpaired(multiplicity)
        args += [ "-c", str ( charge ), "-u", str ( unpaired ) ]
        
        
        #.CPUs
        args += [ "-P", str ( self.parallel ) ]


        #.GFN specify parametrisation of GFN-xTB (default = 2)
        if self.gfn == 3:
            args += [ "--gxtb" ]
        else:
            args += [ "--gfn", str ( self.gfn ) ]

        #.Fermi-smearing
        args += [ "--etemp", str ( self.fermi_temp ) ]

        #.acc accuracy for SCC calculation, lower is better (default = 1.0)
        args += [ "--acc", str ( self.acc ) ]

        #.iterations
        args += [ "--iterations", str ( self.iterations ) ]

        #aditional keys? (a plain string; split into separate arguments)
        if  self.keywords:
            args += self.keywords.split ( )

        #.gradients
        args += [ "--grad" ]

        # adding inputfile and coordinates
        args += [ "--input", state.paths["Input"], state.paths["Coord"] ]

        # Run xtb inside the scratch folder via cwd (no os.chdir of the parent
        # process), sending stdout -- what used to be the shell '> XTBJob.log'
        # redirect -- to the output log, with stderr merged in.
        try:
            with open ( state.paths["Output"], "w" ) as outFile:
                subprocess.run ( args, cwd = directory,
                                 stdout = outFile, stderr = subprocess.STDOUT )
        except Exception:
            # xtb failed to launch/crashed; Energy() detects the missing
            # engrad/output and calls SaveErrorFiles.
            return False
        #' > /dev/null 2>&1'
        try:
            #os.rename(os.path.join(directory,'pcgrad') , os.path.join(directory,'XTBJob.pcgrad'))
            try:
                infile  = open(os.path.join(directory,'pcgrad'), 'r')
                outfile = open(os.path.join(directory,'XTBJob.pcgrad'), 'w')
            
                lines = infile.readlines()
                outfile.write(str(len(lines))+'\n')
                for line in lines:
                    outfile.write(line)
                outfile.close()
                infile.close()
            except:
                #infile  = open(os.path.join(directory,'XTBJob.pc'), 'r')
                
                with open(os.path.join(directory,'XTBJob.pc'), "r") as f:
                    fline = f.readline()
                size = int(fline)
                
                outfile = open(os.path.join(directory,'XTBJob.pcgrad'), 'w')
                outfile.write(str(size)+'\n')
                for line in range(size):
                    outfile.write('0.00000 0.00000 0.00000\n')
                outfile.close()
        except:
            pass
        return True
    
    def Execute_old ( self, state, target ):
        """Execute the xtb job."""
        #print(self._attributable['cpus'])
        #try:
            #outFile = open ( state.paths["Output"], "w" )
            
        backup = os.getcwd()
        directory = os.path.dirname(state.paths["Coord"])
        os.chdir(directory)
        outFile = state.paths["Output"]
        """
        # starting xtb exec
        cmd = self.command
        
        
        # charge and multiplicity
        charge       = target.electronicState.charge
        multiplicity = target.electronicState.multiplicity
        unpaired = self.multiplicity_to_unpaired(multiplicity)
        c_and_u  = ' -c {} -u {} '.format(charge, unpaired)            
        cmd += c_and_u
        
        
        #.CPUs
        cpus = ' -P {} '.format(self.parallel)
        cmd += cpus
        
        
        #.GFN specify parametrisation of GFN-xTB (default = 2)
        cpus = ' --gfn {} '.format(self.gfn)
        cmd += cpus
        
        #.Fermi-smearing
        cmd +=  ' --etemp {} '.format(self.fermi_temp)
        
        #.vip
        #cmd +=  ' --vip '

        
        #.json - write xtbout.json file
        #if  self.json:
        #    cmd +=  ' --json ' 
        
        
        #.acc accuracy for SCC calculation, lower is better (default = 1.0)
        cmd +=  ' --acc {} '.format(self.acc)
        
        #.iterations
        cmd +=  ' --iterations {} '.format(self.iterations)

        
        #.lmo requests localization of orbitals
        #if self.lmo:
        #    cmd +=  ' --lmo '
            
        
        #.fukui : calculates Mulliken partial charges from the neutral, 
        #         positive and negatively charged structure and calculates 
        #         Fukui indices.
        #if  self.vfukui:
        #    cmd +=  ' --vfukui ' 
        
        #aditional keys?
        if  self.keywords:
            #for key in self.keywords:
            cmd += ' ' + self.keywords + ' '
        
        
        #.gradients
        cmd +=  ' --grad '


        # adding inputfile
        cmd +=  ' --input ' + state.paths["Input"] + ' '
        
        # adding coordinates / redicting the outputfile
        #cmd += state.paths["Coord"] #+ ' > ' + state.paths["Output"]
        cmd += state.paths["Coord"] + ' > ' + state.paths["Output"] 
        '''
        /home/fernando/programs/xtb-6.6.1/bin/xtb -c 0 -u 0  -P 1  --gfn 1  
        --etemp 300.0  --acc 1.0  --iterations 300  --grad  
        --input /home/fernando/programs/pDynamo3/scratch/XTBScratch/XTBJob.inp /home/fernando/programs/pDynamo3/scratch/XTBScratch/XTBJob.coord > /home/fernando/programs/pDynamo3/scratch/XTBScratch/XTBJob.log
        '''
        #cmd += state.paths["Coord"] + ' > /dev/null 2>&1'
        
        #print(cmd)
        #subprocess.check_call ( [cmd, state.paths["Coord"]], cwd = directory , stderr = outFile, stdout = outFile )
        os.system(cmd)
        #' > /dev/null 2>&1'
        try:
            #os.rename(os.path.join(directory,'pcgrad') , os.path.join(directory,'XTBJob.pcgrad'))
            infile  = open(os.path.join(directory,'pcgrad'), 'r')
            outfile = open(os.path.join(directory,'XTBJob.pcgrad'), 'w')
        
            lines = infile.readlines()
            outfile.write(str(len(lines))+'\n')
            for line in lines:
                outfile.write(line)
            outfile.close()
            infile.close()
        except:
            pass
        return True
        """
        # Build the command as an ARGUMENT LIST (not a shell string). Safer than
        # os.system (no shell parsing) and lets us run with cwd=directory instead
        # of chdir'ing the whole process -- removing the old cwd-leak bug.
        args = [ self.command ]
        
        
        # charge and multiplicity
        charge       = target.electronicState.charge
        multiplicity = target.electronicState.multiplicity
        unpaired = self.multiplicity_to_unpaired(multiplicity)
        args += [ "-c", str ( charge ), "-u", str ( unpaired ) ]
        
        
        #.CPUs
        args += [ "-P", str ( self.parallel ) ]


        #.GFN specify parametrisation of GFN-xTB (default = 2)
        if self.gfn == 3:
            args += [ "--gxtb" ]
        else:
            args += [ "--gfn", str ( self.gfn ) ]

        #.Fermi-smearing
        args += [ "--etemp", str ( self.fermi_temp ) ]

        #.acc accuracy for SCC calculation, lower is better (default = 1.0)
        args += [ "--acc", str ( self.acc ) ]

        #.iterations
        args += [ "--iterations", str ( self.iterations ) ]

        #aditional keys? (a plain string; split into separate arguments)
        if  self.keywords:
            args += self.keywords.split ( )

        #.gradients
        args += [ "--grad" ]

        # adding inputfile and coordinates
        args += [ "--input", state.paths["Input"], state.paths["Coord"] ]

        # Run xtb inside the scratch folder via cwd (no os.chdir of the parent
        # process), sending stdout -- what used to be the shell '> XTBJob.log'
        # redirect -- to the output log, with stderr merged in.
        try:
            with open ( state.paths["Output"], "w" ) as outFile:
                subprocess.run ( args, cwd = directory,
                                 stdout = outFile, stderr = subprocess.STDOUT )
        except Exception:
            # xtb failed to launch/crashed; Energy() detects the missing
            # engrad/output and calls SaveErrorFiles.
            return False
        #' > /dev/null 2>&1'
        try:
            #os.rename(os.path.join(directory,'pcgrad') , os.path.join(directory,'XTBJob.pcgrad'))
            try:
                infile  = open(os.path.join(directory,'pcgrad'), 'r')
                outfile = open(os.path.join(directory,'XTBJob.pcgrad'), 'w')
            
                lines = infile.readlines()
                outfile.write(str(len(lines))+'\n')
                for line in lines:
                    outfile.write(line)
                outfile.close()
                infile.close()
            except:
                #infile  = open(os.path.join(directory,'XTBJob.pc'), 'r')
                
                with open(os.path.join(directory,'XTBJob.pc'), "r") as f:
                    fline = f.readline()
                size = int(fline)
                
                outfile = open(os.path.join(directory,'XTBJob.pcgrad'), 'w')
                outfile.write(str(size)+'\n')
                for line in range(size):
                    outfile.write('0.00000 0.00000 0.00000\n')
                outfile.close()
        except:
            pass
        return True

    def OrbitalEnergies ( self, target ):
        """Orbital energies and HOMO and LUMO indices."""
        return ( target.scratch.XTBOutputData.get ( "Orbital Energies", None ) ,
                 target.scratch.XTBOutputData.get ( "HOMO"            , -1   ) ,
                 target.scratch.XTBOutputData.get ( "LUMO"            , -1   ) )

    def ReadEngradFile ( self, target, XTBOutputData, gradients3 ):
        """Read an engrad file."""
        # . The energy and gradients are in atomic units.
        state = getattr ( target, self.__class__._stateName )
        try:
            egFile = open ( state.paths["EnGrad"], "r" )
            for i in range ( 7 ): next ( egFile )     # . The number of atoms section and the energy header.
            XTBOutputData["Energy"] = float ( next ( egFile ) )
            for i in range ( 3 ): ( next ( egFile ) ) # . The gradients header.
            for i in range ( len ( state.atomicNumbers ) ):
                for j in range ( 3 ):
                    gradients3[i,j] = float ( ( next ( egFile ) ) )
            egFile.close ( )
            return True
        except:
            return False

    def ReadOutputFile_old (self, target, XTBOutputData):
        """ Function doc """
        state  = getattr ( target, self.__class__._stateName )
        #print(state.paths["Output"])
        atFile = open ( state.paths["Output"], "r" )
        scratch         = { "Is Successful" : False }
        try:
            n = len ( state.atomicNumbers )
            for line in atFile:
                #print(line)
                if line == "Chelpg Charges":
                    data = Array.WithExtent ( n )
                    line = next ( outFile )
                    for i in range ( n ):
                        words   = next ( outFile ).split ( ":", 1 )
                        data[i] = float ( words[-1] )
                    scratch["CHELPG Charges"] = data
                
                elif 'Mulliken/CM5' in line.split():
                    #print(line)
                    data1 = Array.WithExtent ( n )
                    data2 = Array.WithExtent ( n )
                    #line = next ( atFile )
                    for i in range ( n ):
                        #.something like:
                        #['1N', '-0.53592', '-1.03886', '1.367', '4.169', '0.000']
                        words   = next ( atFile ).split ()
                        #print(words)
                        data1[i] = float ( words[1] )
                        data2[i] = float ( words[2] )
                    scratch["Mulliken Charges"] = data1
                    scratch["CM5 Charges"] = data2
                    #print (data1)
                    #print (data2)
                # . Convergence OK if xTB being used (added by Fernando Bachega).
                elif  "convergence criteria satisfied after" in line:
                    words                   = line.split ( )
                    scratch["Cycles"]       = int ( words[5] )
                    scratch["Is Converged"] = True
    
                elif "TOTAL ENERGY" in line:
                    #print(line.split ( ))
                    words                   = line.split ( )
                    scratch["Energy"]       = float( words[3] )
                
                elif "HOMO-LUMO GAP" in line:
                    words                   = line.split ( )
                    scratch["HOMO-LUMO"]    = float( words[3] )
                    
            XTBOutputData.update(scratch)  
            atFile.close ( )
            return True    
        except:
            return False

    def ReadOutputFile ( self, target, XTBOutputData ):
        """Parse the xTB text output into XTBOutputData.

        Modular design: the whole log is read once into memory, then each
        property is extracted by its own independent helper. To stop collecting
        a property, simply comment out its line in the `extractors` list below --
        the helpers do not depend on one another, so disabling one never affects
        the others. Each helper is defensive: it catches its own errors and
        simply skips (leaving its key unset) if its section is absent or its
        format is unexpected, so one malformed section never aborts the rest.
        """
        state = getattr ( target, self.__class__._stateName )
        try:
            with open ( state.paths["Output"], "r" ) as atFile:
                lines = atFile.readlines ( )
        except Exception:
            return False

        n    = len ( state.atomicNumbers )
        data = { "Is Successful" : True }

        # ---- control panel -------------------------------------------------
        # Comment out any line to stop extracting that property. Order is free.
        extractors = [
            self._xtb_extract_energy            ,   # "Energy"            (Eh)
            self._xtb_extract_homo_lumo_gap     ,   # "HOMO-LUMO"         (eV)
            self._xtb_extract_gradient_norm     ,   # "Gradient Norm"     (Eh/a0)
            self._xtb_extract_convergence       ,   # "Cycles", "Is Converged"
            self._xtb_extract_energy_breakdown  ,   # "Energy Terms"      (dict, Eh)
            self._xtb_extract_total_charge      ,   # "Total Charge"      (e)
            self._xtb_extract_partial_charges   ,   # "Partial Charges"   (per atom, any GFN)
            self._xtb_extract_mulliken_cm5      ,   # "Mulliken Charges"/"CM5 Charges" (GFN1 only)
            self._xtb_extract_chelpg_charges    ,   # "CHELPG Charges"    (if present)
            self._xtb_extract_orbitals          ,   # "Orbitals" (+ "HOMO"/"LUMO")
            self._xtb_extract_wiberg_bonds      ,   # "Wiberg Bonds"      (list of (i,j,order))
            self._xtb_extract_dipole            ,   # "Dipole"            (x,y,z,tot Debye)
            self._xtb_extract_metadata          ,   # "XTB Version", "GFN", "Wall Time"
            self._xtb_extract_termination       ,   # "Normal Termination" (bool) + warning
        ]
        # --------------------------------------------------------------------

        for extractor in extractors:
            try:
                extractor ( lines, n, data )
            except Exception:
                # a single failing property must not abort the others
                pass

        XTBOutputData.update ( data )

        # Surface a termination warning to the log so the user notices a run
        # that did not finish normally (the value is set by
        # _xtb_extract_termination). Also mirror it to the pDynamo logFile if
        # available.
        warning = data.get ( "Warning" )
        if warning:
            print ( warning )
            try:
                logFile.Paragraph ( warning )
            except Exception:
                pass

        return True

    # -- individual, independent extractors ---------------------------------
    #    Each fills one (or a few closely-related) key(s) in `data`. They take
    #    (lines, n, data): the full list of log lines, the atom count, and the
    #    output dict to populate. A helper that cannot find its section leaves
    #    its key(s) unset.

    @staticmethod
    def _xtb_find ( lines, needle, start = 0 ):
        """Return the index of the first line containing needle, or -1."""
        for i in range ( start, len ( lines ) ):
            if needle in lines[i]:
                return i
        return -1

    def _xtb_extract_energy ( self, lines, n, data ):
        i = self._xtb_find ( lines, "TOTAL ENERGY" )
        if i >= 0:
            data["Energy"] = float ( lines[i].split ( )[3] )

    def _xtb_extract_homo_lumo_gap ( self, lines, n, data ):
        i = self._xtb_find ( lines, "HOMO-LUMO GAP" )
        if i >= 0:
            data["HOMO-LUMO"] = float ( lines[i].split ( )[3] )

    def _xtb_extract_gradient_norm ( self, lines, n, data ):
        i = self._xtb_find ( lines, "GRADIENT NORM" )
        if i >= 0:
            data["Gradient Norm"] = float ( lines[i].split ( )[3] )

    def _xtb_extract_convergence ( self, lines, n, data ):
        i = self._xtb_find ( lines, "convergence criteria satisfied after" )
        if i >= 0:
            data["Cycles"]       = int ( lines[i].split ( )[5] )
            data["Is Converged"] = True
        elif self._xtb_find ( lines, "convergence criteria cannot be satisfied" ) >= 0:
            data["Is Converged"] = False

    def _xtb_extract_energy_breakdown ( self, lines, n, data ):
        """The SUMMARY block: SCC, dispersion, repulsion, ES/XC terms, etc."""
        start = self._xtb_find ( lines, "SUMMARY" )
        if start < 0:
            return
        terms = {}
        for line in lines[start:start+16]:
            # rows look like ':: SCC energy   -301.79...  Eh ::'
            if "::" not in line:
                continue
            body = line.replace ( "::", " " ).strip ( )
            parts = body.split ( )
            label_words, value = [], None
            for p in parts:
                try:
                    value = float ( p ); break
                except ValueError:
                    label_words.append ( p )
            if value is not None and label_words:
                label = " ".join ( label_words )
                if label.lower ( ) not in ( "summary", ):
                    terms[label] = value
        if terms:
            data["Energy Terms"] = terms

    def _xtb_extract_total_charge ( self, lines, n, data ):
        i = self._xtb_find ( lines, "total charge" )
        if i >= 0:
            for p in lines[i].replace ( "::", " " ).split ( ):
                try:
                    data["Total Charge"] = float ( p ); break
                except ValueError:
                    continue

    def _xtb_extract_partial_charges ( self, lines, n, data ):
        """Per-atom partial charges from the '# Z covCN q C6AA a(0)' table.

        This table is printed for ANY GFN (0/1/2), so it is the robust source of
        atomic charges -- unlike the 'Mulliken/CM5' section which only GFN1
        prints.
        """
        i = self._xtb_find ( lines, "covCN" )
        if i < 0:
            return
        charges = Array.WithExtent ( n )
        got = 0
        for line in lines[i+1:]:
            if line.strip ( ) == "":
                break
            words = line.split ( )
            # expected: idx  Z  sym  covCN  q  C6AA  alpha  -> 7 columns
            if len ( words ) < 7:
                continue
            try:
                idx = int ( words[0] ) - 1
                charges[idx] = float ( words[4] )
                got += 1
            except ( ValueError, IndexError ):
                continue
        if got == n:
            data["Partial Charges"] = charges
            # keep a Mulliken alias so existing AtomicCharges() calls work in
            # every GFN, not only GFN1
            data.setdefault ( "Mulliken Charges", charges )

    def _xtb_extract_mulliken_cm5 ( self, lines, n, data ):
        """The 'Mulliken/CM5' table -- printed by GFN1 only."""
        i = self._xtb_find ( lines, "Mulliken/CM5" )
        if i < 0:
            return
        mull = Array.WithExtent ( n )
        cm5  = Array.WithExtent ( n )
        got = 0
        for k in range ( n ):
            j = i + 1 + k
            if j >= len ( lines ):
                break
            words = lines[j].split ( )
            if len ( words ) < 3:
                break
            try:
                mull[k] = float ( words[1] )
                cm5[k]  = float ( words[2] )
                got += 1
            except ValueError:
                break
        if got == n:
            data["Mulliken Charges"] = mull   # overrides the partial-charge alias
            data["CM5 Charges"]      = cm5

    def _xtb_extract_chelpg_charges ( self, lines, n, data ):
        """CHELPG charges, if the section is present."""
        i = self._xtb_find ( lines, "Chelpg Charges" )   # substring => tolerant of newline/indent
        if i < 0:
            return
        charges = Array.WithExtent ( n )
        got = 0
        for k in range ( n ):
            j = i + 2 + k
            if j >= len ( lines ):
                break
            words = lines[j].split ( ":", 1 )
            try:
                charges[k] = float ( words[-1] )
                got += 1
            except ValueError:
                break
        if got == n:
            data["CHELPG Charges"] = charges

    def _xtb_extract_orbitals ( self, lines, n, data ):
        """Orbital energies/occupations, and the HOMO/LUMO energies (eV)."""
        i = self._xtb_find ( lines, "Occupation" )
        if i < 0:
            return
        orbitals = []
        homo = lumo = None
        for line in lines[i+1:]:
            if "---" in line:
                if orbitals:      # second rule = end of the table
                    break
                continue
            if line.strip ( ) == "":
                break
            words = line.split ( )
            if len ( words ) < 3:
                continue
            tag = ""
            if "(HOMO)" in line: tag = "(HOMO)"
            if "(LUMO)" in line: tag = "(LUMO)"
            numeric = [ w for w in words if w not in ( "(HOMO)", "(LUMO)" ) ]
            try:
                energy_ev = float ( numeric[-1] )
            except ( ValueError, IndexError ):
                continue
            orbitals.append ( energy_ev )
            if tag == "(HOMO)": homo = energy_ev
            if tag == "(LUMO)": lumo = energy_ev
        if orbitals:
            data["Orbitals"] = orbitals
        if homo is not None: data["HOMO"] = homo
        if lumo is not None: data["LUMO"] = lumo

    def _xtb_extract_wiberg_bonds ( self, lines, n, data ):
        """Largest Wiberg bond orders -> list of (i, j, order), 1-based indices."""
        i = self._xtb_find ( lines, "Wiberg" )
        if i < 0:
            return
        # find the table header ('# Z sym total ...'), then skip the dashed rule
        header = self._xtb_find ( lines, "total", i )
        if header < 0:
            return
        start = header + 1
        # skip the '-----' rule line(s) that follow the header
        while start < len ( lines ) and "---" in lines[start]:
            start += 1
        bonds = []
        atom_i = None
        for line in lines[start:]:
            if "---" in line:
                break
            if line.strip ( ) == "":
                continue
            if "--" in line:
                left, right = line.split ( "--", 1 )
                lwords = left.split ( )
                try:
                    atom_i = int ( lwords[0] )
                except ( ValueError, IndexError ):
                    atom_i = None
                rest = right
            else:
                rest = line
            if atom_i is None:
                continue
            rwords = rest.split ( )
            k = 0
            while k + 2 < len ( rwords ):
                try:
                    atom_j = int ( rwords[k] )
                    order  = float ( rwords[k+2] )
                    bonds.append ( ( atom_i, atom_j, order ) )
                except ( ValueError, IndexError ):
                    pass
                k += 3
        if bonds:
            data["Wiberg Bonds"] = bonds

    def _xtb_extract_dipole ( self, lines, n, data ):
        """Molecular dipole: (x, y, z, total) in Debye, from the 'full:' row."""
        i = self._xtb_find ( lines, "molecular dipole" )
        if i < 0:
            return
        for line in lines[i:i+5]:
            if line.strip ( ).startswith ( "full:" ):
                words = line.split ( )
                try:
                    data["Dipole"] = [ float ( w ) for w in words[1:5] ]   # x,y,z,tot
                except ( ValueError, IndexError ):
                    pass
                break

    def _xtb_extract_metadata ( self, lines, n, data ):
        """Version, GFN method and wall-time -- handy for logs/diagnostics."""
        i = self._xtb_find ( lines, "xtb version" )
        if i >= 0:
            words = lines[i].split ( )
            try:
                data["XTB Version"] = words[words.index ( "version" ) + 1]
            except ( ValueError, IndexError ):
                pass
        i = self._xtb_find ( lines, "--gfn" )
        if i >= 0:
            words = lines[i].split ( )
            try:
                data["GFN"] = words[words.index ( "--gfn" ) + 1]
            except ( ValueError, IndexError ):
                pass
        i = self._xtb_find ( lines, "wall-time" )
        if i >= 0:
            data["Wall Time"] = lines[i].split ( ":", 1 )[-1].strip ( )

    def _xtb_extract_termination ( self, lines, n, data ):
        """Check whether xtb terminated normally.

        Sets data["Normal Termination"] to True/False and, when it looks like an
        abnormal run, sets data["Warning"] to a message the caller can surface.

        Robustness note: the exact wording varies by xtb version and run mode.
        A run is treated as FAILED only when a positive error signal is present
        ('abnormal termination of xtb', '[ERROR]', 'ERROR STOP'). It is treated
        as OK when a success signal is present ('normal termination of xtb') OR
        when the run clearly reached the end ('finished run on', which some
        versions print instead of 'normal termination'). This avoids a false
        'abnormal termination' warning on logs (like GFN2 single points) that
        finish cleanly yet never print the literal 'normal termination' line.
        """
        text = "".join ( lines )

        # positive evidence of failure
        failed = ( "abnormal termination of xtb" in text ) \
                 or ( "[ERROR]" in text ) \
                 or ( "ERROR STOP" in text )

        # positive evidence of success
        succeeded = ( "normal termination of xtb" in text ) \
                    or ( "finished run on" in text )

        if failed:
            data["Normal Termination"] = False
            data["Warning"] = "warning: xTB abnormal termination"
        elif succeeded:
            data["Normal Termination"] = True
        else:
            # no clear signal either way -- flag conservatively so the user can
            # check, but say it is unconfirmed rather than definitely abnormal
            data["Normal Termination"] = False
            data["Warning"] = "warning: xTB termination could not be confirmed"





  
    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCModelXTB, self ).SummaryItems ( )
        if self.keywords is not None:
            n = len ( self.keywords )
            items.append ( ( "Keywords", "{:s}".format ( "/".join ( self.keywords[0:min(2,n)] ) ) ) )
        return items

    # . Coordinates are written in atomic units (keywords must contain "Bohrs").
    def WriteInputFile ( self, target, doGradients, doQCMM, coordinates3 ):
        """Write an input file."""
        state  = getattr ( target, self.__class__._stateName )
        
        
        coordFile = open ( state.paths["Coord"], "w" )
        coordFile.write ( "$coord\n")
        for ( i, n ) in enumerate ( state.atomicNumbers ):
            coordFile.write ( "{:20.10f}{:20.10f}{:20.10f} {:<12s}\n".format ( 
                                                                            coordinates3[i,0] ,
                                                                            coordinates3[i,1] ,
                                                                            coordinates3[i,2] ,
                                                                            PeriodicTable.Symbol ( n ))
                                                                            )

        coordFile.write ( "$end\n")
        coordFile.close ( )
        
        inFile = open ( state.paths["Input"], "w" )
        
        if doQCMM     : 
            #inFile.write ( '%pointcharges "' + state.paths["PC"] + '"\n' )
            #inFile.write ( '%pointcharges "' + state.paths["PC"] + '"\n' )
            
            inFile.write ( '$embedding\n')
            inFile.write ( '    interface=orca\n'.format(state.paths["PC"]))
            inFile.write ( '    input={}\n'.format(state.paths["PC"]))
            inFile.write ( '$end\n')
        
        
        #inFile = open ( state.paths["Input"], "w" )
        
        #inFile.write ( "#\n" )
        #inFile.write ( "# XTB Job.\n" )
        #inFile.write ( "#\n" )
        
        #if doGradients: mode = "ENGRAD"
        #else          : mode = "ENERGY"
        #inFile.write ( "! " + mode + " BOHRS " + " ".join ( self.keywords ) + "\n" )
        #inFile.write ( "* xyz {:d} {:d}\n".format ( target.electronicState.charge, target.electronicState.multiplicity ) )
        #for ( i, n ) in enumerate ( state.atomicNumbers ):
        #    inFile.write ( "{:<12s}{:20.10f}{:20.10f}{:20.10f}\n".format ( PeriodicTable.Symbol ( n ) ,
        #                                                                   coordinates3[i,0]          ,
        #                                                                   coordinates3[i,1]          ,
        #                                                                   coordinates3[i,2]        ) )
        #inFile.write ( "*\n" )
        #inFile.close ( )

    @property
    def command ( self ):
        """Get the command to execute the program."""
        command = self.__dict__.get ( "_command", None )
        if command is None:
            command = os.getenv ( _XTBCommand )
            # . Command must point to an executable file.
            if  ( command is None ) or not ( os.path.isfile ( command ) and os.access ( command, os.X_OK ) ):
                raise NotInstalledError ( "XTB executable not found." )
            else:
                self.__dict__["_command"] = command
        return command

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

