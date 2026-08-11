"""The XTB QC model."""

import glob, math, os, os.path, subprocess, re, tempfile, shutil

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

# . Scratch base-directory environment variable.
_XTBScratchEnv = "PDYNAMO3_SCRATCH"

# . Scratch directory.
#   [PORTABILITY] Built from PDYNAMO3_SCRATCH, but that variable may be unset on
#   this machine (e.g. a system prepared elsewhere is opened here). Guard against
#   os.getenv returning None so merely importing this module never fails; the
#   real, machine-local resolution happens in _resolve_scratch() at run time.
_scratch_base = os.getenv ( _XTBScratchEnv )
if _scratch_base:
    _XTBScratch = os.path.join ( _scratch_base, "XTBScratch" )
else:
    _XTBScratch = os.path.join ( tempfile.gettempdir ( ), "XTBScratch" )

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
                # The scratch folder is now always a unique random directory
                # (see DeterminePaths/_resolve_scratch), so it is safe to remove
                # it whole. rmtree (unlike the old os.rmdir) also clears files
                # that are NOT prefixed by the job name -- notably 'xtbrestart',
                # whose leftover from a previous system was the source of the
                # cross-system "Index out of range" error.
                scratch = self.paths.get ( "Scratch", None )
                if scratch is not None and os.path.isdir ( scratch ):
                    shutil.rmtree ( scratch, ignore_errors = True )
                else:
                    # fallback: remove only this job's files
                    jobFiles = glob.glob ( os.path.join ( self.paths["Glob"] + "*" ) )
                    for jobFile in jobFiles: os.remove ( jobFile )
            except:
                pass

    @staticmethod
    def _resolve_scratch ( scratch ):
        """Resolve a usable scratch base directory on THIS machine.

        Portability + isolation logic (mirrors the executable fallback):
          1. if the scratch directory assigned to the object exists (or its
             parent exists so it can be created), use it -- this is the folder
             configured where the system was prepared;
          2. otherwise fall back to this machine's PDYNAMO3_SCRATCH;
          3. otherwise fall back to the system temporary directory.

        In every case a UNIQUE random subfolder is then created underneath, so
        two processes / systems can never accidentally share a scratch folder
        (which is what let a residual 'xtbrestart' from another system be read
        and produce an "Index out of range" error).

        Returns the absolute path of the freshly created, unique scratch folder.
        """
        def _usable ( base ):
            if not base:
                return False
            # usable if it already exists, or its parent exists (so we can mkdir)
            if os.path.isdir ( base ):
                return True
            parent = os.path.dirname ( os.path.normpath ( base ) )
            return bool ( parent ) and os.path.isdir ( parent )

        # 1) the path stored on the object (from where the system was prepared)
        base = scratch if _usable ( scratch ) else None

        # 2) fall back to this machine's PDYNAMO3_SCRATCH
        if base is None:
            env_base = os.getenv ( _XTBScratchEnv )
            if env_base:
                env_base = os.path.join ( env_base, "XTBScratch" )
                if _usable ( env_base ):
                    base = env_base

        # 3) last resort: the system temp directory
        if base is None:
            base = os.path.join ( tempfile.gettempdir ( ), "XTBScratch" )

        # make sure the base exists, then create a unique random subfolder
        os.makedirs ( base, exist_ok = True )
        unique = os.path.join ( base, RandomString ( ) )
        os.makedirs ( unique, exist_ok = True )
        return unique

    def DeterminePaths ( self, scratch, deleteJobFiles = True, randomJob = False, randomScratch = False ):
        """Determine the paths needed by an XTB job."""
        paths = {}
        if randomJob: job = RandomString ( )
        else:         job = _DefaultJobName
        # [PORTABILITY + ISOLATION] Always resolve the scratch on THIS machine
        # and run inside a unique random subfolder. _resolve_scratch() handles
        # the case where the assigned path came from another computer and does
        # not exist here, and the random subfolder prevents cross-process/
        # cross-system contamination (e.g. a stale 'xtbrestart'). Because the
        # folder is always unique, it is always safe to remove afterwards, so we
        # record it under "Scratch" for DeleteJobFiles.
        scratch          = self._resolve_scratch ( scratch )
        paths["Scratch"] = scratch
        if not os.path.exists ( scratch ): os.makedirs ( scratch, exist_ok = True )
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
        # Expose the actual (unique, random) scratch folder resolved for this
        # run on the model itself, so code outside pDynamo can find where the
        # job files really are. Needed because the files now live in a random
        # subfolder under self.scratch, not directly in self.scratch -- e.g.
        # EasyHybrid's backup_xtb_files() copies the xTB log from here.
        self.__dict__["_activeScratch"] = state.paths.get ( "Scratch", self.scratch )
        return state

    @property
    def activeScratch ( self ):
        """The actual scratch folder used by the most recent build/run.

        This is the unique random subfolder created by DeterminePaths, i.e.
        where the xTB job files (log, engrad, ...) for the current run really
        live. Falls back to the configured scratch base if no run has been
        built yet.
        """
        return self.__dict__.get ( "_activeScratch", self.scratch )

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

    def Execute(self, state, target):
        """Execute the xTB job."""

        directory = os.path.dirname(state.paths["Coord"])

        charge = target.electronicState.charge
        multiplicity = target.electronicState.multiplicity
        unpaired = self.multiplicity_to_unpaired(multiplicity)

        cmd = [
            self.command,
            "-c", str(charge),
            "-u", str(unpaired),
            "-P", str(self.parallel),
        ]

        # GFN parametrization.
        if self.gfn == 3:
            cmd.append("--gxtb")
        else:
            cmd.extend(["--gfn", str(self.gfn)])

        # Electronic temperature / Fermi smearing.
        cmd.extend(["--etemp", str(self.fermi_temp)])

        # SCC convergence / iterations.
        cmd.extend(["--acc", str(self.acc)])
        cmd.extend(["--iterations", str(self.iterations)])

        # Additional xTB keywords.
        if self.keywords:
            cmd.extend(self.keywords.split())

        # Request gradients.
        cmd.append("--grad")

        # xTB input file.
        cmd.extend(["--input", state.paths["Input"]])

        # Coordinate file.
        cmd.append(state.paths["Coord"])

        # Execute xTB.
        with open(state.paths["Output"], "w") as output:
            result = subprocess.run(
                cmd,
                cwd=directory,
                stdout=output,
                stderr=subprocess.STDOUT,
                check=False
            )

        if result.returncode != 0:
            return False

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


        ## Process point-charge gradients.
        #pcgrad = os.path.join(directory, "pcgrad")
        #output_pcgrad = os.path.join(directory, "XTBJob.pcgrad")
        #
        #if os.path.exists(pcgrad):
        #    with open(pcgrad, "r") as infile:
        #        lines = infile.readlines()
        #
        #    with open(output_pcgrad, "w") as outfile:
        #        outfile.write(str(len(lines)) + "\n")
        #        outfile.writelines(lines)
        #
        #else:
        #    # No pcgrad produced by xTB.
        #    pcfile = os.path.join(directory, "XTBJob.pc")
        #
        #    if not os.path.exists(pcfile):
        #        return False
        #
        #    with open(pcfile, "r") as f:
        #        size = int(f.readline())
        #
        #    with open(output_pcgrad, "w") as outfile:
        #        outfile.write(str(size) + "\n")
        #        for _ in range(size):
        #            outfile.write("0.00000 0.00000 0.00000\n")
        #
        #return True
    
    def Execute_old ( self, state, target ):
        """Execute the xtb job."""
        #print(self._attributable['cpus'])
        #try:
            #outFile = open ( state.paths["Output"], "w" )
            
        backup = os.getcwd()
        directory = os.path.dirname(state.paths["Coord"])
        os.chdir(directory)
        outFile = state.paths["Output"]
        
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
        if self.gfn == 3:
            cpus = ' --gxtb '.format(self.gfn)
        else:
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
        
        #print('target',target)
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
        ]
        # --------------------------------------------------------------------

        for extractor in extractors:
            try:
                extractor ( lines, n, data )
            except Exception:
                # a single failing property must not abort the others
                pass

        XTBOutputData.update ( data )
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
        #print("Energy", data["Energy"])

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
        
    @staticmethod
    def _is_valid_executable ( command ):
        """True if 'command' is a non-empty path to an existing executable file."""
        return ( command is not None ) and os.path.isfile ( command ) \
               and os.access ( command, os.X_OK )

    @property
    def command ( self ):
        """Get the command to execute the program.

        Resolution order (makes an exported system portable across machines):
          1. the path stored on this object (self.__dict__["_command"]), if it
             still points to a valid executable on THIS machine -- this is the
             path configured where the system was prepared;
          2. otherwise, the PDYNAMO3_XTBCOMMAND environment variable of the
             CURRENT machine -- so a system prepared elsewhere (with a different
             xTB path) still runs here;
          3. if neither is valid, raise NotInstalledError.

        The old behaviour trusted a cached '_command' blindly: a system prepared
        on machine A carried A's absolute xTB path, and on machine B it kept
        trying that non-existent path instead of falling back to B's env var.
        Validating the path before use fixes that.
        """
        # 1) path stored on the object (from where the system was prepared)
        stored = self.__dict__.get ( "_command", None )
        if self._is_valid_executable ( stored ):
            return stored

        # 2) fall back to this machine's environment variable
        env_command = os.getenv ( _XTBCommand )
        if self._is_valid_executable ( env_command ):
            # remember the valid path so we don't re-check every call
            self.__dict__["_command"] = env_command
            return env_command

        # 3) neither worked
        raise NotInstalledError (
            "XTB executable not found. Checked the path stored with the system "
            "({}) and the {} environment variable ({}). Set {} to the xtb "
            "executable on this machine.".format (
                stored, _XTBCommand, env_command, _XTBCommand ) )

    def SetCommand ( self, path, validate = True ):
        """Redefine the xTB executable path used by this model.

        Stores 'path' as the executable to run (self.__dict__["_command"]),
        taking precedence over the PDYNAMO3_XTBCOMMAND environment variable on
        the next run.

        By default the path is validated (must be an existing executable file);
        an invalid path raises NotInstalledError so the caller finds out
        immediately instead of at calculation time. Pass validate=False to store
        the path unconditionally (e.g. when configuring a machine where the file
        is not present yet). Pass path=None to clear the stored path and fall
        back to the environment variable again.
        """
        if path is None:
            self.__dict__.pop ( "_command", None )
            return
        if validate and not self._is_valid_executable ( path ):
            raise NotInstalledError (
                "Not a valid xTB executable: {}".format ( path ) )
        self.__dict__["_command"] = path

    def SetScratch ( self, path, validate = True ):
        """Redefine the scratch base directory used by this model.

        Stores 'path' as the scratch base (self.scratch). It is resolved on the
        actual machine at run time by _resolve_scratch (which also creates a
        unique random subfolder under it), so a value that does not exist yet is
        fine -- it will be created, or fall back to PDYNAMO3_SCRATCH / the system
        temp directory if it turns out not to be usable here.

        By default the path is lightly validated: it must be creatable, i.e. it
        already exists or its parent directory exists. Pass validate=False to
        store it unconditionally. Pass path=None to reset to the default scratch.
        """
        if path is None:
            self.scratch = _XTBScratch
            return
        if validate:
            usable = os.path.isdir ( path ) or \
                     os.path.isdir ( os.path.dirname ( os.path.normpath ( path ) ) )
            if not usable:
                raise QCModelError (
                    "Scratch directory is not creatable (neither it nor its "
                    "parent exists): {}".format ( path ) )
        self.scratch = path

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

