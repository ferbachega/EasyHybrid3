#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Lembrar de colocar uma header nesse arquivo

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import glob, math, os, os.path, sys

from datetime import date
import time

import numpy as np

from pprint import pprint

#VISMOL_HOME = os.environ.get('VISMOL_HOME')

#path fo the core python files on your machine
#sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/pDynamoMethods") )
#sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/gui"))

#from LogFile import LogFileReader

#from gEngine.vismol_object import EVismolObject

#from vismol.model.atom import Atom
from vismol.model.residue import Residue
from vismol.model.chain import Chain
from vismol.core.vismol_object import VismolObject

from logging import getLogger
logger = getLogger(__name__)

#---------------------------------------
from pBabel                    import*                                     
from pCore                     import*  
#---------------------------------------
from pMolecule                 import*                              
from pMolecule.MMModel         import*
from pMolecule.NBModel         import*                                     
from pMolecule.QCModel         import*
#---------------------------------------
from pScientific               import*                                     
from pScientific.Arrays        import*                                     
from pScientific.Geometry3     import*                                     
from pScientific.RandomNumbers import*                                     
from pScientific.Statistics    import*
from pScientific.Symmetry      import*
#---------------------------------------                              
from pSimulation               import*
#---------------------------------------

'''
#import our core lib
from SimulationsPreset import Simulation 
#---------------------------------------
from vModel import VismolObject
from vModel.MolecularProperties import ATOM_TYPES_BY_ATOMICNUMBER
from vModel.MolecularProperties import COLOR_PALETTE

'''

'''
from easyhybrid.gui import *
from easyhybrid.gui.PES_analisys_window  import  PotentialEnergyAnalysisWindow 
from easyhybrid.gui.PES_analisys_window  import  parse_2D_scan_logfile 
'''

class pDynamoSession:
    """ Class doc """
    
    def __init__ (self, vm_session = None):
        """ Class initialiser """
        self.vm_session  = vm_session
        self.name        = 'p_session'
        
        self.nbModel_default         = NBModelCutOff.WithDefaults ( )
        self.fixed_color             = [0.5, 0.5, 0.5]
        self.pdynamo_distance_safety = 0.5
        
        '''self.active_id is the attribute that tells which 
        system is active for calculations in pdynamo 
        (always an integer value)'''
        self.active_id       = 0
        
        
        '''Now we can open more than one pdynamo system. 
        Each new system loaded into memory is stored in 
        the form of a dictionary, which has an int as 
        an access key.'''
        self.psystem =  {
                         0:None
                        }
        
        #self.systems_list = []
        self.counter      = 0
        self.color_palette_counter = 0
        #self.sel_name_counter = 0


    def load_a_new_pDynamo_system_from_dict (self, input_files = {}, system_type = 0, name = None, tag = None):
        """ Function doc """

        #psystem = self.generate_pSystem_dictionary()
        
        system = None 
        if system_type == 0:
            system              = ImportSystem       ( input_files['amber_prmtop'] )
            system.coordinates3 = ImportCoordinates3 ( input_files['coordinates'] )
            self.define_NBModel(_type = 1, system = system)                      
        
        elif system_type == 1:
            parameters          = CHARMMParameterFileReader.PathsToParameters (input_files['charmm_par'])
            system              = ImportSystem       ( input_files['charmm_psf'] , isXPLOR = True, parameters = parameters )
            system.coordinates3 = ImportCoordinates3 ( input_files['coordinates'] )
            self.define_NBModel(_type = 1, system = system)        
        
        elif system_type == 2:
            mmModel        = MMModelOPLS.WithParameterSet ( input_files['opls_folder'] )       
            system         = ImportSystem ( input_files['coordinates'])
            system.DefineMMModel ( mmModel )
            self.define_NBModel(_type = 1, system = system)          
        
        elif system_type == 3 or system_type == 4 :
            system = ImportSystem (input_files['coordinates'])

        else:
            pass
        system.Summary()
        
        
        system.e_input_files = input_files
        
        
        if tag:
            system.e_tag = tag
        else:
            system.e_tag ='MolSys'
        
        if name:
            system.label = name
        else:
            pass
        
        
        self._add_attributes_to_pSystem (system, working_folder = None)
        
        self.append_system_to_pdynamo_session ( system)


    def append_system_to_pdynamo_session (self, system):
        """ Function doc """

        
        try:
            system.e_charges_backup           = list(system.AtomicCharges()).copy()
        except:
            system.e_charges_backup           = []
        
        #psystem['system']                  =  system
        #psystem['name']                    =  name
        #psystem['tag']                     =  tag
        #psystem['color_palette']           =  COLOR_PALETTE[self.color_palette_counter]
        
        system.e_id                      =  self.counter
        self.psystem[system.e_id]        =  system 
        
        #print('color_palette', self.color_palette_counter)
        #self.systems_list.append(psystem)
        
        self.active_id   = self.counter  
        self.counter    += 1

        
        # obtaing the color palette
        '''
        if self.color_palette_counter >= len(COLOR_PALETTE)-1:
            self.color_palette_counter = 0
        else:
            self.color_palette_counter += 1
        '''
        

        self._build_vobject_from_pDynamo_system ( system = system ) 
        
        '''
        self.systems[self.active_id]['step_counter'] += 1
        if self.vm_session.main_session.selection_list_window.visible:
            self.vm_session.main_session.refresh_system_liststore()
            self.vm_session.main_session.selection_list_window.update_window(system_names = True, coordinates = False,  selections = False)
        '''


    def _add_attributes_to_pSystem (self, system, working_folder = None ):
        """ Function doc """

        system.e_id                       = 0               # access number (same as the access key in the self.psystem dictionary)
        #system.e_tag                      = tag             # 15 length string

        system.e_date                     = date.today()    # Time     
        system.e_color_palette            = None            # will be replaced by a dict

        #system.e_vobject                  = None            # Vismol object associated with the system -> is the object that will 

        system.e_active                   = False          
        system.e_bonds                    = None           
        system.e_sequence                 = None           


        system.e_qc_table                 = None           
        system.e_fixed_table              = []             
        system.e_charges_backup           = []             
        system.e_selections               = {}             

        system.e_vobjects                 = {}                  
        system.e_logfile_data             = {}             # <--- vobject_id : [data0, data1, data2], ...]  obs: each "data" is a dictionary 
        system.e_working_folder           = working_folder 
        system.e_step_counter             = 0              


    def _get_sequence_from_pDynamo_system (self, system = None):
        """ Function doc """
        
        if system:
            sequence = getattr ( system, "sequence", None )
            
            if sequence is None: 
                sequence = Sequence.FromAtoms (system.atoms, componentLabel = "UNK.1" )
            else:
                pass
        else:
            sequence = None
        
        return sequence


    def _get_atom_info_from_pdynamo_atom_obj (self, sequence = None, atom = None):
        """
        To extract information from atom objects, 
        belonging to pdynamo, and to  organize it as a list
        that will be delivered later to build the 
        vismolObj
        
        """
        entityLabel = atom.parent.parent.label
        useSegmentEntityLabels = False
        if useSegmentEntityLabels:
            chainID = ""
            segID   = entityLabel[0:4]
        else:
            chainID = entityLabel[0:1]
            segID   = ""

        
        if sequence:
            resName, resSeq, iCode = sequence.ParseLabel ( atom.parent.label, fields = 3 )
        else:
            sequence = None
        
        ##print(atom.index, atom.label,resName, resSeq, iCode,chainID, segID,  atom.atomicNumber, atom.connections)#, xyz[0], xyz[1], xyz[2] )
        
        index        = atom.index
        at_name      = atom.label
        at_resi      = int(resSeq)
        at_resn      = resName
        at_ch        = chainID
        try:
            at_symbol    = ATOM_TYPES_BY_ATOMICNUMBER[atom.atomicNumber] # at.get_symbol(at_name)
        except:
            at_symbol = "O"
        
        at_occup     = 0.0
        at_bfactor   = 0.0
        at_charge    = 0.0
        atom         = {
              'index'      : index      , 
              'name'       : at_name    , 
              'resi'       : at_resi    , 
              'resn'       : at_resn    , 
              'chain'      : at_ch      , 
              'symbol'     : at_symbol  , 
              'occupancy'  : at_occup   , 
              'bfactor'    : at_bfactor , 
              'charge'     : at_charge   
              }
        
        return atom


    def _get_qc_table_from_pDynamo_system(self, system):
        system.e_qc_table = list(system.qcState.pureQCAtoms)


    def _build_vobject_from_pDynamo_system (self                                          , 
                                            system                    = None              , 
                                            vm_object_name            = 'a_new_vismol_obj',
                                            e_id                      = None              ,
                                            is_vobject_active         = True              ,
                                            autocenter                = True              ,
                                            #refresh_qc_and_fixed     = True              ,
                                            #add_vobject_to_vm_session = True             ,
                                            frames                    = None
                                           ):
        '''
        atoms.append({
                      'index'      : index      , 
                      'name'       : at_name    , 
                      'resi'       : at_resi    , 
                      'resn'       : at_resn    , 
                      'chain'      : at_ch      , 
                      'symbol'     : at_symbol  , 
                      'occupancy'  : at_occup   , 
                      'bfactor'    : at_bfactor , 
                      'charge'     : at_charge   
                      })
        
        '''
        sequence = self._get_sequence_from_pDynamo_system(system)
        #print (sequence)


        atoms = []     
        frame = []
        
        
        atom_qtty = len(system.atoms.items)
        coords    = np.empty([1, atom_qtty, 3], dtype=np.float32)
        j = 0
        for atom in system.atoms.items:
            xyz = system.coordinates3[atom.index]
            x = np.float32(xyz[0])
            y = np.float32(xyz[1])
            z = np.float32(xyz[2])
            coords[0,j,:] = x, y, z
            #frame.append(xyz[0])
            #frame.append(xyz[1])
            #frame.append(xyz[2])
            atoms.append(self._get_atom_info_from_pdynamo_atom_obj(sequence =  sequence, atom   = atom))
            j += 1
            
        #frame = np.array(frame, dtype=np.float32)
      
        #print(atoms)
        #if frames is None:
        #    frames = [frame]
        #else:
        #    pass
        
        
        #atom_qtty = len(system.atoms.items)
        #
        #coords = np.empty([len(frames), atom_qtty, 3], dtype=np.float32)
        #for atom in system.atoms.items:
        ##for frame in frames:
        #    j = 0
        #    lines = frame.strip().split("\n")
        #    for line in lines:
        #        if line[:4] == "ATOM" or line[:6] == "HETATM":
        #            x = np.float32(line[30:38])
        #            y = np.float32(line[38:46])
        #            z = np.float32(line[46:54])
        #            coords[i,j,:] = x, y, z
        #            j += 1
        #    i += 1
        #print("_get_frame_coords(frames, atom_qtty)")
        #print(coords)
        #return coords        
                
        
        
        
        
        
        vm_object = VismolObject(self.vm_session, 
                                 len(self.vm_session.vm_objects_dic), 
                                 name=vm_object_name)
                                 
        vm_object.set_model_matrix(self.vm_session.vm_glcore.model_mat)
        
        #mudar isso depois
        unique_id = len(self.vm_session.atom_dic_id)
        
        initial = time.time()
        
        atom_id = 0
        pprint(atoms)
        for _atom in atoms:
            if _atom["chain"] not in vm_object.chains.keys():
                vm_object.chains[_atom["chain"]] = Chain(vm_object, name=_atom["chain"])
            _chain = vm_object.chains[_atom["chain"]]
            
            if _atom["resi"] not in _chain.residues.keys():
                _r = Residue(vm_object, name=_atom["resn"], index=_atom["resi"], chain=_chain)
                vm_object.residues[_atom["resi"]] = _r
                _chain.residues[_atom["resi"]] = _r
            _residue = _chain.residues[_atom["resi"]]
            
            print(_atom)
            #atom = Atom()
            atom = Atom(vismol_object = vm_object,#
                        name          =_atom["name"] , 
                        index         =_atom["index"],
                        residue       =_residue      , 
                        chain         =_chain, 
                        atom_id       =atom_id,
                        occupancy     =_atom["occupancy"], 
                        bfactor       =_atom["bfactor"],
                        charge        =_atom["charge"])
                        

                        
                        
            atom.unique_id = unique_id
            atom._generate_atom_unique_color_id()
            _residue.atoms[atom_id] = atom
            vm_object.atoms[atom_id] = atom
            atom_id += 1
            unique_id += 1
        
        logger.debug("Time used to build the tree: {:>8.5f} secs".format(time.time() - initial))
        
        
        
        #coords = np.empty([0, atom_id, 3], dtype=np.float32)
        #for f in frames:
        #    coords = np.vstack((coords, f))
            
        vm_object.frames = coords
        vm_object.mass_center = np.mean(vm_object.frames[0], axis=0)
        vm_object.find_bonded_and_nonbonded_atoms()
        
        vm_object._generate_color_vectors(self.vm_session.atom_id_counter)
        vm_object.active = True
        self.vm_session._add_vismol_object(vm_object, show_molecule=True)
        return vm_object





    def define_a_new_QCModel (self, system = None, parameters = None):
        """ Function doc """
        
        '''Here we have to reload the mmModel original charges.
        This is postulated because multiple associations of QC 
        regions can distort the charge distribution of some residues. 
        (because charge rescale)'''
        
        electronicState = ElectronicState.WithOptions ( charge = parameters['charge'], multiplicity = parameters['multiplicity'])
        system.electronicState = electronicState
        
        converger = DIISSCFConverger.WithOptions( energyTolerance   = parameters['energyTolerance'] ,
                                                  densityTolerance  = parameters['densityTolerance'] ,
                                                  maximumIterations = parameters['maximumIterations'] )


        qcModel         = QCModelMNDO.WithOptions ( hamiltonian = parameters['method'], converger=converger )
        #_QCmodel = QCModelMNDO.WithOptions( hamiltonian = _method, converger=converger )

        
        if system.e_qc_table:
            
            '''This function rescales the remaining charges in the MM part. The 
            sum of the charges in the MM region must be an integer value!'''
            self.check_charge_fragmentation()
            '''----------------------------------------------------------------'''
            
            system.DefineQCModel (qcModel, qcSelection = Selection.FromIterable ( system.e_qc_table) )
            #self.refresh_qc_and_fixed_representations()# static = False )
            
            if system.mmModel:
                system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
            else:
                pass
            #self.add_a_new_item_to_selection_list (system_id = self.active_id, 
            #                                         indexes = self.systems[self.active_id]['qc_table'], 
            #                                         name    = 'QC atoms')
            
            #print(self.systems[self.active_id]['selections']["QC atoms"])
        else:
            system.DefineQCModel (qcModel)
            
            #self.refresh_qc_and_fixed_representations()
            
            
            #self.add_a_new_item_to_selection_list (system_id = self.active_id, 
            #                                         indexes = range(0, self.systems[self.active_id]['system'].atoms), 
            #                                         name    = 'QC atoms')
            #    

    def check_charge_fragmentation(self, correction = True):
        """ Function doc """
        #self.systems[self.active_id]['system_original_charges']
        mm_residue_table = {}
        qc_residue_table = self.systems[self.active_id]['qc_residue_table']
        system           = self.systems[self.active_id]['system']
        
        print('\n\n\Sum of total charges(MM)', sum(system.mmState.charges))
        '''----------------------------------------------------------------'''
        '''Restoring the original charges before rescheduling a new region.'''
        original_charges = self.systems[self.active_id]['system_original_charges'].copy()
        
        for index, charge in enumerate(original_charges):
            system.mmState.charges[index]   = original_charges[index]
        '''----------------------------------------------------------------'''
        print('\n\n\Sum of total charges(MM)', sum(system.mmState.charges))

        qc_charge        = 0.0
        
        if system.mmModel is None:
            return None 
        
        '''Here we are going to arrange the atoms that are not in the QC part, 
        but are in the same residues as those atoms within the QC part.'''  

        self._check_ref_vobject_in_pdynamo_system()
        
        for res in self.systems[self.active_id]['vobject'].residues:
            
            if res.resi in qc_residue_table.keys():
                
                mm_residue_table[res.resi] = []
                
                for atom in res.atoms:
                    index_v = atom.index-1
                    index_p = system.atoms.items[index_v].index
                    index_p = system.atoms.items[index_v].label
                    #print( system.mmState.charges)
                    charge  = system.mmState.charges[index_v]
                    resn    = res.resn 
                    atom.charge = system.mmState.charges[index_v]
                    
                    if index_v in qc_residue_table[res.resi]:
                        qc_charge += atom.charge
                        pass
                        #print (resn, res.resi, index_v, index_p, charge, True )
                    
                    else:
                        #print (resn, res.resi, index_v, index_p, charge, False)
                        mm_residue_table[res.resi].append(index_v)
                
                
                
                #print(atom.index, atom.atomicNumber, system.mmState.charges[idx],self.systems[self.active_id]['vobject'].atoms[idx].resn )
            
        #print('mm_residue_table',mm_residue_table)
        '''Here we are going to do a rescaling of the charges of atoms of 
        the MM part but that the residues do not add up to an entire charge.''' 
        
        
        
        for resi in mm_residue_table.keys():
            
            total = 0.0
            for index in mm_residue_table[resi]:
                pcharge = system.mmState.charges[index]
                total += pcharge
            
            rounded  = float(round(total))
            diff     = rounded - total
            size     = len(mm_residue_table[resi])
            
            if size > 0:
                fraction = diff/size
                print('residue: ', resi, 'charge fraction = ',fraction)
            
                for index in mm_residue_table[resi]:
                    system.mmState.charges[index] += fraction
                    #total += pcharge
            else:
                pass
        print('\n\n\Sum of total charges(MM) after rescaling', sum(system.mmState.charges))
        print('\n\n\Sum of total charges(MM) original', sum(self.systems[self.active_id]['system_original_charges']))
        
        print('QC charge from selected atoms: ',round(qc_charge) )
        #for atom in self.systems[self.active_id]['vobject'].atoms:
        #    print( atom.index, atom.name, atom.charge)
        #print('Total charge: ', sum(system.mmState.charges))










    def build_vobject_from_pDynamo_system (self                                         , 
                                           name                     = 'a_new_vismol_obj',
                                           #system_id                = None              ,
                                           #vobject_active           = True              ,
                                           #autocenter               = True              ,
                                           #refresh_qc_and_fixed     = True              ,
                                           #add_vobject_to_vm_session = True              ,
                                           #frames                   = None
                                           ):
        
        """ Function doc """
        print('build_vobject_from_pDynamo_system')
        
        #'''
        if system_id is not None:
            pass
        else:
            system_id = self.active_id
        
        
        #self.get_bonds_from_pDynamo_system(safety = self.pdynamo_distance_safety, system_id = system_id)
        self.get_sequence_from_pDynamo_system(system_id = system_id)

        atoms = []     
        frame = []
        
        for atom in self.systems[system_id]['system'].atoms.items:
            xyz = self.get_atom_coords_from_pdynamo_system (atom   = atom, system  = self.systems[system_id]['system'])
            frame.append(xyz[0])
            frame.append(xyz[1])
            frame.append(xyz[2])
            atoms.append(self.get_atom_info_from_pdynamo_atom_obj(atom   = atom, system_id = system_id))
        frame = np.array(frame, dtype=np.float32)
      
        
        if frames is None:
            frames = [frame]
        else:
            pass


        name  = os.path.basename(name)

        vobject  = VismolObject.VismolObject(name                           = name                                          ,    
                                             atoms                          = atoms                                         ,    
                                             vm_session                     = self.vm_session                            ,    
                                             #bonds_pair_of_indexes          = self.systems[system_id]['bonds']         ,    
                                             trajectory                     = frames                                       ,  
                                             color_palette                  = self.systems[system_id]['color_palette'] ,
                                             auto_find_bonded_and_nonbonded = True               )

        vobject.easyhybrid_system_id = self.systems[system_id]['id']
        vobject.set_model_matrix(self.vm_session.glwidget.vm_widget.model_mat)
        vobject.active = vobject_active
        vobject._get_center_of_mass(frame = 0)
        if self.systems[system_id]['system'].qcModel:
            sum_x = 0.0 
            sum_y = 0.0 
            sum_z = 0.0
            
            self.systems[system_id]['qc_table'] = list(self.systems[system_id]['system'].qcState.pureQCAtoms)
            total = 0
            
            for atom_index in self.systems[system_id]['qc_table']:
                atom = vobject.atoms[atom_index]
                
                coord = atom.coords (frame = 0)
                sum_x += coord[0]
                sum_y += coord[1]
                sum_z += coord[2]
                total+=1
                
                    
            center = np.array([sum_x / total,
                               sum_y / total, 
                               sum_z / total])
            
        else:
            center = vobject.mass_center

        self.systems[system_id]['vobject'] = vobject
        if add_vobject_to_vm_session:
            self.vm_session.add_vobject_to_vismol_session (pdynamo_session  = self,
                                                                    #rep              = True, 
                                                                    vobject    = vobject, 
                                                                    autocenter       = autocenter)
        
        if refresh_qc_and_fixed:
            self.refresh_qc_and_fixed_representations(system_id = system_id, vobject = vobject) 

        self.vm_session.glwidget.vm_widget.center_on_coordinates(vobject, center)
        self.vm_session.main_session.update_gui_widgets()
        return vobject
        #'''





















    def get_system_name (self, system_id = None):
        """ Function doc """
        if system_id:
            return self.psystem[system_id]['name']
        else:
            return self.psystem[self.active_id]['name']




    def define_NBModel (self, _type = 1 , parameters =  None, system = None):
        """ Function doc """
        
        if _type == 0:
            self.nbModel = NBModelFull.WithDefaults ( )
        
        elif _type == 1:
            self.nbModel = NBModelCutOff.WithDefaults ( )
        
        elif _type == 2:
            self.nbModel = NBModelORCA.WithDefaults ( )
        
        elif _type == 3:
            self.nbModel = NBModelDFTB.WithDefaults ( )
        
        
        try:
            if system:
                system.DefineNBModel ( self.nbModel )
                system.Summary ( )
            else:
                self.systems[self.active_id]['system'].DefineNBModel ( self.nbModel )
                self.systems[self.active_id]['system'].Summary ( )
            return True
        
        except:
            print('failed to bind nbModel')
            return False
        










import numpy as np
from vismol.model.molecular_properties import ATOM_TYPES


class Atom:
    """ Class doc """
    
    def __init__(self, vismol_object=None, name="Xx", index=None, residue=None,
                 chain=None, pos=None, symbol=None, atom_id=None, color=None,
                 vdw_rad=None, cov_rad=None, ball_rad=None,
                 occupancy=0.0, bfactor=0.0, charge=0.0, bonds_indexes=None):
        """ Class initializer """
        self.vm_object = vismol_object
        self.name = name
        self.index = index   # - Remember that the "index" attribute refers to the numbering of atoms (it is not a zero base, it starts at 1 for the first atom)
        self.residue = residue
        self.chain = chain
        self.pos = pos     # - coordinates of the first frame
        self.unique_id = None
        
        # self.symbol = self._get_symbol(self.name) if (symbol is None) else symbol
        if symbol is None:
            self.symbol = self._get_symbol()
        else:
            self.symbol = symbol
        self.atom_id = atom_id
        
        if color is None:
            self.color = self._init_color()
        else:
            self.color = color
        
        if vdw_rad is None:
            self.vdw_rad = self._init_vdw_rad()
        else:
            self.vdw_rad = vdw_rad
        
        if cov_rad is None:
            self.cov_rad = self._init_cov_rad()
        else:
            self.cov_rad = cov_rad
        
        if ball_rad is None:
            self.ball_rad = self._init_ball_rad()
        else:
            self.ball_rad = ball_rad
        
        self.color_id = None
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.charge = charge
        if bonds_indexes is None:
            self.bonds_indexes = []
        else:
            self.bonds_indexes = bonds_indexes
        
        self.selected       = False
        self.lines          = True
        self.dots           = False
        self.nonbonded      = False
        self.impostor       = False
        self.ribbons        = False
        self.ball_and_stick = False
        self.sticks         = False
        self.spheres        = False
        self.dash           = False
        self.surface        = False
        self.bonds          = []
        self.isfree         = True
    
    def _get_symbol(self):
        """ Function doc """
        name = self.name.strip()
        if name == "":
            return ""
        
        _n = name
        for char in name:
            if char.isnumeric():
                _n = _n.replace(char, "")
        name = _n
        
        if len(name) >= 3:
            name = name[:2]
        
        # This can't happen since before all numbers have been converted to strings
        # if len(name) == 2:
        #     if name[1].isnumeric():
        #         symbol = name[0]
        
        # The capitalization of hte name can solve all the next ifs
        # name = name.lower().capitalize()
        if name in ATOM_TYPES.keys():
            return name
        else:
            if name[0] == "H":
                if name[1] == "g":
                    symbol =  "Hg"
                elif name[1] =="e":
                    symbol = "He"
                else:
                    symbol =  "H"
            
            elif name[0] == "C":
                if name[1] == "a":
                    symbol = "Ca"
                elif name[1] =="l":
                    symbol = "Cl"
                elif name[1] =="d":
                    symbol = "Cd"
                elif name[1] =="u":
                    symbol = "Cu"
                else:
                    symbol = "C"
            
            elif name[0] == "N":
                if name[1] == "i":
                    symbol = "Ni"
                elif name[1] == "a":
                    symbol = "Na"
                elif name[1] == "e":
                    symbol = "Ne"
                elif name[1] == "b":
                    symbol = "Nb"
                else:
                    symbol = "N"
            
            elif name[0] == "O":
                if name[1] == "s":
                    symbol = "Os"
                else:
                    symbol = "O"
            
            elif name[0] == "S":
                if name[1] == "I":
                    symbol = "Si"
                elif name[1] == "e":
                    symbol = "Se"
                else:
                    symbol = "S"
            
            elif name[0] == "P":
                if name[1] == "d":
                    symbol = "Pd"
                elif  name[1] == "b":
                    symbol = "Pb"
                elif  name[1] == "o":
                    symbol = "Po"
                else:
                    symbol = "P" 
            
            elif name[0] == "M":
                if name[1] == "n":
                    symbol = "Mn"
                elif name[1] == "N":
                    symbol = "Mn"
                elif name[1] == "o":
                    symbol = "Mo"
                elif name[1] == "G":
                    symbol = "Mg"
                else:
                    symbol = "X"
            else:
                symbol = "X"
        return symbol
    
    def _init_color(self):
        """ Return the color of an atom in RGB. Note that the returned
            value is in scale of 0 to 1, but you can change this in the
            index. If the atomname does not match any of the names
            given, it returns the default dummy value of atom X.
        """
        try:
            color = self.vm_object.color_palette[self.name]
        except KeyError:
            # print(self.symbol, "Atom")
            color = self.vm_object.color_palette[self.symbol]
        return np.array(color, dtype=np.float32)
    
    def _generate_atom_unique_color_id(self):
        """ Function doc """
        r = (self.unique_id & 0x000000FF) >> 0
        g = (self.unique_id & 0x0000FF00) >> 8
        b = (self.unique_id & 0x00FF0000) >> 16
        # pickedID = r + g * 256 + b * 256*256
        self.color_id = np.array([r/255.0, g/255.0, b/255.0], dtype=np.float32)
    
    def _init_vdw_rad(self):
        """ Function doc """
        try:
            vdw = ATOM_TYPES[self.name][6]
        except KeyError:
            vdw = ATOM_TYPES[self.symbol][6]
        return vdw
    
    def _init_cov_rad(self):
        """ Function doc """
        try:
            cov = ATOM_TYPES[self.name][5]
        except KeyError:
            cov = ATOM_TYPES[self.symbol][5]
        return cov
    
    def _init_ball_rad(self):
        """ Function doc """
        try:
            ball = ATOM_TYPES[self.name][6]
        except KeyError:
            ball = ATOM_TYPES[self.symbol][6]
        return ball
    
    def coords(self, frame=None):
        """ Function doc """
        if frame is None:
            frame  = self.vm_object.vm_session.frame
        # x = self.vm_object.frames[frame][(self.index-1)*3  ]
        # y = self.vm_object.frames[frame][(self.index-1)*3+1]
        # z = self.vm_object.frames[frame][(self.index-1)*3+2]
        # return np.array([x, y, z])
        return self.vm_object.frames[frame, self.atom_id]
    
    def get_grid_position(self, gridsize=3, frame=None):
        """ Function doc """
        coords = self.coords(frame)
        gridpos = (int(coords[0]/gridsize), int(coords[1]/gridsize), int(coords[2]/gridsize))
        return gridpos
    
    def get_cov_rad(self):
        """ Function doc """
        return self.cov_rad
    
    def init_color_rgb(self, name):
        """ Return the color of an atom in RGB. Note that the returned
            value is in scale of 0 to 1, but you can change this in the
            index. If the atomname does not match any of the names
            given, it returns the default dummy value of atom X.
        """
        try:
            color = color =self.vm_object.color_palette[name]
            #color = ATOM_TYPES[name][1]
        except KeyError:
            if name[0] == "H":# or name in self.hydrogen:
                #color = ATOM_TYPES["H"][1]
                color = self.vm_object.color_palette["H"]
            
            elif name[0] == "C":
                #color = ATOM_TYPES["C"][1]
                color = self.vm_object.color_palette["C"]
            
            elif name[0] == "O":
                #color = ATOM_TYPES["O"][1]
                color = self.vm_object.color_palette["O"]
            
            elif name[0] == "N":
                #color = ATOM_TYPES["N"][1]
                color = self.vm_object.color_palette["N"]
                
            elif name[0] == "S":
                #color = ATOM_TYPES["S"][1]
                color = self.vm_object.color_palette["S"]
            else:
                #color = ATOM_TYPES["X"][1]
                color = self.vm_object.color_palette["X"]
                
        color = [int(color[0]*250), int(color[1]*250), int(color[2]*250)]
        return color
    
    def init_radius(self, name):
        """
        """
        try:
            rad = ATOM_TYPES[name][6]/5.0
        except KeyError:
            if name[0] == "H" or name in self.hydrogen:
                rad = ATOM_TYPES["H"][6]/5.0
            elif name[0] == "C":
                rad = ATOM_TYPES["C"][6]/5.0
            elif name[0] == "O":
                rad = ATOM_TYPES["O"][6]/5.0
            elif name[0] == "N":
                rad = ATOM_TYPES["N"][6]/5.0
            elif name[0] == "S":
                rad = ATOM_TYPES["S"][6]/5.0
            else:
                rad = 0.30
        return rad






