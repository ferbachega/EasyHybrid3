"""The NWChem QC model.

Drives NWChem as an external QC engine for pDynamo/EasyHybrid, with QM/MM
support via NWChem's Bq module (external point charges loaded from a file).

STATUS / WHAT IS VALIDATED
--------------------------
* Structure, scratch/command portability, subprocess execution: modelled on the
  (working) QCModelXTB and safe.
* Input generation (WriteInputFile): follows the official NWChem documentation
  (geometry/basis/dft/charge/bq/task blocks). Confident but should be checked
  against a real run.
* Output parser (ReadOutputFile / ReadGradient): based on documented NWChem
  output patterns ("Total DFT energy", "ENERGY GRADIENTS"). The exact header
  offsets MUST be confirmed against a REAL NWChem output file -- these are the
  lines most likely to need small adjustments. They are marked [VALIDATE] below.
* QM/MM point-charge file: NWChem 'bq load <file> format 1 2 3 4 units bohr'
  reads x y z q. The exact format of the point-charge file pDynamo produces
  must be confirmed and, if needed, rewritten in WriteInputFile. Marked
  [VALIDATE:QMMM].
"""

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
# . Default error prefix.
_DefaultErrorPrefix = "error_"

# . Default job name.
_DefaultJobName = "NWChemJob"

# . Command environment variable.
_NWChemCommand = "PDYNAMO3_NWCHEMCOMMAND"

# . Scratch base-directory environment variable.
_NWChemScratchEnv = "PDYNAMO3_SCRATCH"

# . Scratch directory (import-safe: PDYNAMO3_SCRATCH may be unset on this box).
_scratch_base = os.getenv ( _NWChemScratchEnv )
if _scratch_base:
    _NWChemScratch = os.path.join ( _scratch_base, "NWChemScratch" )
else:
    _NWChemScratch = os.path.join ( tempfile.gettempdir ( ), "NWChemScratch" )

#===================================================================================================================================
# . State class.
#===================================================================================================================================
class QCModelNWChemState ( QCModelState ):
    """A QC model state."""

    _attributable = dict ( QCModelState._attributable )
    _attributable.update ( { "deleteJobFiles" : False ,
                             "paths"          : None  } )

    def __del__ ( self ):
        """Deallocation."""
        self.DeleteJobFiles ( )

    def DeleteJobFiles ( self ):
        """Delete job files."""
        if self.deleteJobFiles:
            try:
                scratch = self.paths.get ( "Scratch", None )
                if scratch is not None and os.path.isdir ( scratch ):
                    shutil.rmtree ( scratch, ignore_errors = True )
                else:
                    jobFiles = glob.glob ( os.path.join ( self.paths["Glob"] + "*" ) )
                    for jobFile in jobFiles: os.remove ( jobFile )
            except:
                pass

    def DeterminePaths ( self, scratch, deleteJobFiles = False, randomJob = False, randomScratch = False ):
        """Determine the paths needed by an NWChem job."""
        paths = {}
        if randomJob: job = RandomString ( )
        else:         job = _DefaultJobName
        if randomScratch:
            scratch          = os.path.join ( scratch, RandomString ( ) )
            paths["Scratch"] = scratch  # . Only set (and later removed) if random.
        if not os.path.exists ( scratch ):
            os.makedirs ( scratch, exist_ok = True )
        jobRoot       = os.path.join ( scratch, job )
        paths["Glob"] = jobRoot
        for ( key, ext ) in ( ( "Input"  , "nw"  ) ,   # NWChem input deck
                              ( "Output" , "out" ) ,   # NWChem text output
                              ( "PC"     , "pc"  ) ):   # external point charges
            paths[key] = "{:s}.{:s}".format ( jobRoot, ext )
        self.deleteJobFiles = deleteJobFiles
        self.paths          = paths

    def SaveErrorFiles ( self, message ):
        """Save the input and output files for inspection if there is an error."""
        for key in ( "Input", "Output" ):
            path = self.paths.get ( key, None )
            if path and os.path.exists ( path ):
                ( head, tail ) = os.path.split ( path )
                os.rename ( path, os.path.join ( head, _DefaultErrorPrefix + tail ) )
        ( head, tail ) = os.path.split ( self.paths["Glob"] )
        raise QCModelError ( message + "\nCheck the files \"{:s}*\".".format ( os.path.join ( head, _DefaultErrorPrefix + tail ) ) )

#===================================================================================================================================
# . Model class.
#===================================================================================================================================
class QCModelNWChem ( QCModel ):
    """The NWChem QC model."""

    _attributable = dict ( QCModel._attributable )
    _classLabel   = "NWChem QC Model"
    _stateName    = "nwchemState"
    _stateObject  = QCModelNWChemState
    _summarizable = dict ( QCModel._summarizable )
    _attributable.update ( { "deleteJobFiles" : False        ,
                             "randomJob"      : False        ,
                             "randomScratch"  : False        ,
                             "scratch"        : _NWChemScratch,
                             # NWChem-specific options:
                             "functional"     : "b3lyp"      ,   # dft xc
                             "basis"          : "6-31G*"     ,   # basis library
                             "method"         : "dft"        ,   # dft | scf | mp2 ...
                             "charge"         : 0            ,
                             "keywords"       : None         ,   # extra dft-block lines (list)
                             "memory"         : None         } ) # e.g. "2000 mb"
    _summarizable.update ( { "functional"     : "Functional"       ,
                             "basis"          : "Basis Set"        ,
                             "method"         : "Method"           } )

    def AtomicCharges ( self, target, chargeModel = ChargeModel.Mulliken ):
        """Atomic charges."""
        source = getattr ( target.scratch, "nwchemOutputData", {} )
        return source.get ( "Mulliken Charges", None )

    def AtomicSpins ( self, target, chargeModel = "Mulliken" ):
        """Atomic spins."""
        return None

    def BondOrders ( self, target, chargeModel = None ):
        """Bond orders."""
        source = getattr ( target.scratch, "nwchemOutputData", {} )
        return source.get ( "Bond Orders", None )

    def BuildModel ( self, target, qcSelection = None ):
        """Build the model."""
        state = super ( QCModelNWChem, self ).BuildModel ( target, qcSelection = qcSelection )
        state.DeterminePaths ( self.scratch                         ,
                               deleteJobFiles = self.deleteJobFiles ,
                               randomJob      = self.randomJob      ,
                               randomScratch  = self.randomScratch  )
        return state

    def DipoleMoment ( self, target, center = None ):
        """Dipole moment."""
        source = getattr ( target.scratch, "nwchemOutputData", {} )
        return source.get ( "Dipole", None )

    def Energy ( self, target ):
        """Calculate the quantum chemical energy."""
        doGradients      = target.scratch.doGradients
        nwchemOutputData = {}
        state            = getattr ( target, self.__class__._stateName )
        target.scratch.nwchemOutputData = nwchemOutputData

        self.WriteInputFile ( target, doGradients, ( target.nbModel is not None ), target.scratch.qcCoordinates3AU )

        isOK = self.Execute ( state, target )
        if not isOK: state.SaveErrorFiles ( "Error executing program." )

        if doGradients:
            isOK = self.ReadGradient ( target, nwchemOutputData, target.scratch.qcGradients3AU )
            if not isOK: state.SaveErrorFiles ( "Error reading gradient." )

        isOK = self.ReadOutputFile ( target, nwchemOutputData )
        if not isOK: state.SaveErrorFiles ( "Error reading output file." )

        target.scratch.energyTerms["NWChem QC"] = ( nwchemOutputData["Energy"] * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )

    def multiplicity_to_mult ( self, multiplicity = 1 ):
        """NWChem's 'mult' is the spin multiplicity directly (2S+1)."""
        return multiplicity

    def Execute ( self, state, target ):
        """Execute the NWChem job via subprocess (no shell, no chdir of parent)."""
        directory = os.path.dirname ( state.paths["Input"] )
        if not os.path.isdir ( directory ):
            try:    os.makedirs ( directory, exist_ok = True )
            except Exception: pass

        # NWChem is invoked as:  nwchem <input.nw>   with stdout -> output file.
        # (Parallel runs would prefix mpirun -np N; left out here -- see note.)
        args = [ self.command, state.paths["Input"] ]
        try:
            with open ( state.paths["Output"], "w" ) as outFile:
                subprocess.run ( args, cwd = directory,
                                 stdout = outFile, stderr = subprocess.STDOUT )
        except Exception:
            return False
        return True

    def OrbitalEnergies ( self, target ):
        """Orbital energies and HOMO and LUMO indices."""
        source = getattr ( target.scratch, "nwchemOutputData", {} )
        return ( source.get ( "Orbital Energies", None ) ,
                 source.get ( "HOMO", -1 ) ,
                 source.get ( "LUMO", -1 ) )

    #-------------------------------------------------------------------------------
    # . Gradient reading.
    #-------------------------------------------------------------------------------
    def ReadGradient ( self, target, nwchemOutputData, gradients3 ):
        """Read the energy gradient from the NWChem output.

        [VALIDATE] NWChem prints a block like:

                          DFT ENERGY GRADIENTS

            atom               coordinates                        gradient
                         x          y          z           x          y          z
           1 O       0.000000   0.000000   0.000000    0.000000   0.000000  -0.001234
           ...

        Columns 5,6,7 are the gradient (atomic units, Hartree/Bohr). The header
        offset (how many lines between the 'ENERGY GRADIENTS' title and the first
        atom row) and the exact title string should be confirmed on a real run.
        """
        state = getattr ( target, self.__class__._stateName )
        try:
            with open ( state.paths["Output"], "r" ) as f:
                lines = f.readlines ( )
        except Exception:
            return False

        n = len ( state.atomicNumbers )

        # locate the LAST 'ENERGY GRADIENTS' block (the converged one)
        idx = -1
        for i, line in enumerate ( lines ):
            if "ENERGY GRADIENTS" in line:
                idx = i
        if idx < 0:
            return False

        # [VALIDATE] Instead of assuming a fixed header size, scan forward from
        # the title for the first line that looks like an atom row (integer
        # index + 8 columns) and read n consecutive such rows. This is robust to
        # the exact number of header lines, which varies by module/version.
        got = 0
        j = idx + 1
        while j < len ( lines ) and got < n:
            words = lines[j].split ( )
            j += 1
            if len ( words ) < 8:
                if got > 0:
                    break     # left the table after having started
                continue
            # an atom row starts with an integer index
            try:
                int ( words[0] )
            except ValueError:
                if got > 0:
                    break
                continue
            try:
                gradients3[got,0] = float ( words[-3] )
                gradients3[got,1] = float ( words[-2] )
                gradients3[got,2] = float ( words[-1] )
                got += 1
            except ValueError:
                if got > 0:
                    break
        return got == n

    #-------------------------------------------------------------------------------
    # . Modular output parser (same "control panel" idea as QCModelXTB).
    #-------------------------------------------------------------------------------
    def ReadOutputFile ( self, target, nwchemOutputData ):
        """Parse the NWChem text output into nwchemOutputData.

        Comment a line out of the `extractors` list to stop collecting that
        property. Each helper is defensive and independent.
        """
        state = getattr ( target, self.__class__._stateName )
        try:
            with open ( state.paths["Output"], "r" ) as f:
                lines = f.readlines ( )
        except Exception:
            return False

        n    = len ( state.atomicNumbers )
        data = { "Is Successful" : True }

        extractors = [
            self._nw_extract_energy        ,   # "Energy" (Hartree)
            self._nw_extract_convergence   ,   # "Is Converged"
            self._nw_extract_mulliken      ,   # "Mulliken Charges"
            self._nw_extract_dipole        ,   # "Dipole"
            self._nw_extract_homo_lumo     ,   # "HOMO"/"LUMO"
            self._nw_extract_termination   ,   # "Normal Termination" + warning
        ]

        for extractor in extractors:
            try:
                extractor ( lines, n, data )
            except Exception:
                pass

        nwchemOutputData.update ( data )

        warning = data.get ( "Warning" )
        if warning:
            print ( warning )
            try:    logFile.Paragraph ( warning )
            except Exception: pass

        # Energy is required downstream; fail if it was not found.
        return ( "Energy" in data )

    @staticmethod
    def _nw_find ( lines, needle, start = 0 ):
        for i in range ( start, len ( lines ) ):
            if needle in lines[i]:
                return i
        return -1

    @staticmethod
    def _nw_find_last ( lines, needle ):
        idx = -1
        for i, line in enumerate ( lines ):
            if needle in line:
                idx = i
        return idx

    def _nw_extract_energy ( self, lines, n, data ):
        """Total energy. [VALIDATE] NWChem DFT prints 'Total DFT energy ='.

        Other methods use different labels ('Total SCF energy', 'Total MP2
        energy', ...); we try the DFT one first, then a couple of fallbacks.
        Take the LAST occurrence (converged value).
        """
        labels = [ "Total DFT energy", "Total SCF energy",
                   "Total MP2 energy", "Total CCSD energy" ]
        for label in labels:
            i = self._nw_find_last ( lines, label )
            if i >= 0:
                # line: '         Total DFT energy =      -76.012345678'
                try:
                    data["Energy"] = float ( lines[i].split ( "=" )[-1].split ( )[0] )
                    return
                except ( ValueError, IndexError ):
                    continue

    def _nw_extract_convergence ( self, lines, n, data ):
        """SCF/DFT convergence flag. [VALIDATE] wording varies by module."""
        if self._nw_find ( lines, "calculation not converged" ) >= 0 or \
           self._nw_find ( lines, "did not converge" ) >= 0:
            data["Is Converged"] = False
        elif self._nw_find ( lines, "Total DFT energy" ) >= 0 or \
             self._nw_find ( lines, "Total SCF energy" ) >= 0:
            data["Is Converged"] = True

    def _nw_extract_mulliken ( self, lines, n, data ):
        """Mulliken charges. [VALIDATE] block layout must be confirmed.

        NWChem prints a 'Mulliken analysis of the total density' section; the
        per-atom charge column and offset depend on the module/version, so this
        is intentionally conservative and may need adjusting on a real file.
        """
        i = self._nw_find ( lines, "Mulliken analysis of the total density" )
        if i < 0:
            return
        charges = Array.WithExtent ( n )
        got = 0
        for line in lines[i+1:]:
            words = line.split ( )
            if not words:
                continue
            # rows typically look like: '  1 O    8   8.34  ... charge'
            # [VALIDATE] pick the charge column against a real file
            try:
                # heuristic: an atom row starts with an integer index
                int ( words[0] )
            except ValueError:
                if got > 0:
                    break     # left the table
                continue
            try:
                charges[got] = float ( words[-1] )
                got += 1
            except ( ValueError, IndexError ):
                pass
            if got >= n:
                break
        if got == n:
            data["Mulliken Charges"] = charges

    def _nw_extract_dipole ( self, lines, n, data ):
        """Dipole moment vector. [VALIDATE] label/units to confirm."""
        i = self._nw_find_last ( lines, "Dipole moment" )
        if i < 0:
            return
        # often the components follow on nearby lines as 'DMX DMY DMZ'
        for line in lines[i:i+6]:
            if "DMX" in line or "Dipole moment" in line:
                nums = re.findall ( r'[-+]?\d+\.\d+', line )
                if len ( nums ) >= 3:
                    data["Dipole"] = [ float ( x ) for x in nums[:3] ]
                    return

    def _nw_extract_homo_lumo ( self, lines, n, data ):
        """HOMO/LUMO from orbital energies. [VALIDATE] NWChem 'Vector' listing.

        NWChem prints molecular orbital 'Vector' lines with 'Occ=' and 'E='.
        The HOMO is the last with non-zero occupation; LUMO the first empty one.
        """
        homo = lumo = None
        prev_e = None
        for line in lines:
            if "Occ=" in line and "E=" in line:
                try:
                    occ = float ( line.split ( "Occ=" )[1].split ( )[0].replace ( "D", "E" ) )
                    e   = float ( line.split ( "E=" )[1].split ( )[0].replace ( "D", "E" ) )
                except ( ValueError, IndexError ):
                    continue
                if occ > 0.5:
                    homo = e
                elif lumo is None:
                    lumo = e
        if homo is not None: data["HOMO"] = homo
        if lumo is not None: data["LUMO"] = lumo

    def _nw_extract_termination ( self, lines, n, data ):
        """Check for normal termination.

        NWChem prints a CITATION / 'Total times' block and the line
        'Task  times' at successful completion; failures print 'error' /
        'NWChem input module' aborts. Success signals used (any of):
        'Total times cpu', 'CITATION'. Failure: 'For further details see the...'
        after an error, or 'ERROR'. [VALIDATE] against real successful/failed
        runs.
        """
        text = "".join ( lines )
        failed    = ( "This error has not been assigned" in text ) or \
                    ( "  ERROR  " in text ) or \
                    ( "input module" in text and "error" in text.lower ( ) )
        succeeded = ( "Total times  cpu" in text ) or \
                    ( "Task  times  cpu" in text ) or \
                    ( "CITATION" in text )
        if failed and not succeeded:
            data["Normal Termination"] = False
            data["Warning"] = "warning: NWChem abnormal termination"
        elif succeeded:
            data["Normal Termination"] = True
        else:
            data["Normal Termination"] = False
            data["Warning"] = "warning: NWChem termination could not be confirmed"

    #-------------------------------------------------------------------------------
    # . Input generation.
    #-------------------------------------------------------------------------------
    def SummaryItems ( self ):
        """Summary items."""
        items = [ ( self._classLabel, True ) ,
                  ( "Method",     "{:s}".format ( self.method ) ) ,
                  ( "Functional", "{:s}".format ( str ( self.functional ) ) ) ,
                  ( "Basis Set",  "{:s}".format ( str ( self.basis ) ) ) ]
        return items

    def WriteInputFile ( self, target, doGradients, doQCMM, coordinates3 ):
        """Write an NWChem input deck.

        Layout (units are Bohr, matching pDynamo's qcCoordinates3AU):

            start <job>
            [memory <mem>]
            charge <q>
            geometry units bohr nocenter noautosym noautoz
              <sym> <x> <y> <z>
              ...
            end
            basis
              * library <basis>
            end
            dft
              xc <functional>
              mult <multiplicity>
              [<extra keyword lines>]
            end
            [bq load <pc file> format 1 2 3 4 units bohr]     # QM/MM
            task <method> <gradient|energy>
        """
        state = getattr ( target, self.__class__._stateName )

        charge       = target.electronicState.charge
        multiplicity = target.electronicState.multiplicity
        mult         = self.multiplicity_to_mult ( multiplicity )

        job = _DefaultJobName

        inFile = open ( state.paths["Input"], "w" )
        inFile.write ( "start {:s}\n".format ( job ) )
        inFile.write ( "title \"Generated by EasyHybrid\"\n" )
        if self.memory:
            inFile.write ( "memory {:s}\n".format ( str ( self.memory ) ) )
        inFile.write ( "charge {:d}\n".format ( charge ) )

        # geometry (Bohr; no reorientation so gradients map back to pDynamo atoms)
        inFile.write ( "geometry units bohr nocenter noautosym noautoz\n" )
        for ( i, num ) in enumerate ( state.atomicNumbers ):
            inFile.write ( "  {:<3s} {:20.10f} {:20.10f} {:20.10f}\n".format (
                PeriodicTable.Symbol ( num ) ,
                coordinates3[i,0] ,
                coordinates3[i,1] ,
                coordinates3[i,2] ) )
        inFile.write ( "end\n" )

        # basis
        inFile.write ( "basis\n" )
        inFile.write ( "  * library {:s}\n".format ( str ( self.basis ) ) )
        inFile.write ( "end\n" )

        # method block (DFT shown; SCF/MP2 would differ)
        if self.method.lower ( ) == "dft":
            inFile.write ( "dft\n" )
            inFile.write ( "  xc {:s}\n".format ( str ( self.functional ) ) )
            inFile.write ( "  mult {:d}\n".format ( mult ) )
            if self.keywords:
                for line in self.keywords:
                    inFile.write ( "  {:s}\n".format ( line ) )
            inFile.write ( "end\n" )

        # QM/MM external point charges via Bq
        if doQCMM:
            # [VALIDATE:QMMM] The pc file must contain 'x y z q' rows in Bohr.
            # If pDynamo writes them in a different layout/units, convert here.
            inFile.write ( "bq load {:s} format 1 2 3 4 units bohr\n".format ( state.paths["PC"] ) )

        # task
        mode = "gradient" if doGradients else "energy"
        inFile.write ( "task {:s} {:s}\n".format ( self.method.lower ( ), mode ) )
        inFile.close ( )

    @property
    def command ( self ):
        """Get the NWChem executable, with a portable fallback chain."""
        stored = self.__dict__.get ( "_command", None )
        if self._is_valid_executable ( stored ):
            return stored
        env_command = os.getenv ( _NWChemCommand )
        if self._is_valid_executable ( env_command ):
            self.__dict__["_command"] = env_command
            return env_command
        raise NotInstalledError (
            "NWChem executable not found. Checked the path stored with the "
            "system ({}) and the {} environment variable ({}). Set {} to the "
            "nwchem executable on this machine.".format (
                stored, _NWChemCommand, env_command, _NWChemCommand ) )

    @staticmethod
    def _is_valid_executable ( command ):
        """True if command is a path to an existing executable file."""
        return ( command is not None ) and os.path.isfile ( command ) \
               and os.access ( command, os.X_OK )

    def SetCommand ( self, path, validate = True ):
        """Redefine the NWChem executable path (see QCModelXTB.SetCommand)."""
        if path is None:
            self.__dict__.pop ( "_command", None )
            return
        if validate and not self._is_valid_executable ( path ):
            raise NotInstalledError ( "Not a valid NWChem executable: {}".format ( path ) )
        self.__dict__["_command"] = path

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
