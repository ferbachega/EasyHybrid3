#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  easyhybrid_pDynamo_selection.py
#  
#  Copyright 2022 Fernando <fernando@winter>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
import gc
import os

#import math
from gi.repository import Gtk, Pango
from gi.repository import Gdk
from gi.repository import GdkPixbuf
#from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import FolderChooserButton

from gui.windows.setup.windows_and_dialogs import TextWindow
from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox


from pprint import pprint
import numpy as np

class ReimagingTrajectoryWindow:

    def __init__(self, main = None , system_liststore = None):
        """ Class initialiser """
        self.main                = main
        self.home                = main.home
        self.p_session           = self.main.p_session
        self.vm_session          = self.main.vm_session
        self.Visible             = False  
        
        self.system_liststore      = system_liststore
        self.coordinates_liststore = Gtk.ListStore(str, int)
        
        self.parameters = []
        
        self.types = {
                     0: 'Protein'           ,
                     1: 'N-Ca-C'            ,
                     2: 'C-alpha'           ,
                     3: 'Current Selection' ,
                     }
        

    def OpenWindow (self, vobject = None):
        """ Function doc """
        if self.Visible  ==  False:

            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.main.home,'src/gui/windows/analysis/reimaging_trajectory.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_title('Reimaging Trajectory')
            #self.window.set_keep_above(True)            
            self.window.set_default_size(400, 100)
            self.btn_reimaging = self.builder.get_object('btn_reimaging')
             
            self.coordinates_liststore = Gtk.ListStore(str, int, int)
            self.box_coordinates = self.builder.get_object('box_coordinates')
            #self.combobox_coordinates = CoordinatesComboBox(self.main.vobject_liststore_dict[self.main.p_session.active_id])
            self.coordinates_combobox = CoordinatesComboBox(self.main.vobject_liststore_dict[self.p_session.active_id])
            self.box_coordinates.pack_start  (self.coordinates_combobox, False, False, 0)
            

            self.box_system       = self.builder.get_object('box_system')
            self.combobox_systems = SystemComboBox(self.main, self.coordinates_combobox)
            self.box_system.pack_start  (self.combobox_systems, False, False, 0)

            self.combobox_systems.set_active(0)
            self.coordinates_combobox.set_active(0)
            
            self.method_store = Gtk.ListStore(str)
            methods = [
                       'by geometric center'       ,
                       #'by center of mass'         , 
                       ]
            
            for method in methods:
                self.method_store.append([method])
            self.cb_align_by = self.builder.get_object('cb_align_by')
            self.cb_align_by.set_model(self.method_store)
            
            renderer_text = Gtk.CellRendererText()
            self.cb_align_by.pack_start(renderer_text, True)
            self.cb_align_by.add_attribute(renderer_text, "text", 0)
            self.cb_align_by.set_active(0)
            

            self.btn_reimaging.connect("clicked", self.on_btn_reimaging)
            self.combobox_systems.connect("changed", self.on_combobox_systems_changed)
            
            self.coordinates_combobox.connect("changed", self.on_combobox_coord_changed)
            
            self.window.show_all()           
            self.Visible  = True
            
            
    def on_combobox_systems_changed (self, widget):
        """ Function doc """
        cb_id = widget.get_system_id()
        
        if cb_id is not None:
            
            self.update_window (coordinates = True)
            
            key =  self.coordinates_combobox.get_vobject_id()
            #_, key = self.coordinates_liststore[cb_id]
            
            self.VObj = self.vm_session.vm_objects_dic[key]

            
    def update_window (self, system_names = False, coordinates = True,  selections = False ):
        """ Function doc """

        if self.Visible:
            
            _id = self.combobox_systems.get_active()
            if _id == -1:
                '''_id = -1 means no item inside the combobox'''
                return None
            else:    
                _, system_id, pixbuf = self.system_liststore[_id]
            
            if system_names:
                self.refresh_system_liststore ()
                self.combobox_systems.set_active(_id)
            
            if coordinates:
                self.refresh_coordinates_liststore ()
    

    def on_combobox_coord_changed (self, widget):
        """ Function doc """
        cb_id = widget.get_system_id()
        print('cb',cb_id )

    def refresh_coordinates_liststore(self, system_id = None):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
        #print(2313, system_id,self.main.vobject_liststore_dict )
        self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
        self.coordinates_combobox.set_active_vobject(-1)
            

    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
        #self.on_btn_clear(None, None)

    def on_btn_reimaging (self, widget, event = None):
        """
        Reimage molecules within the periodic boundary conditions (PBC).

        This function checks whether molecules are outside the simulation cell 
        and repositions them back into the central unit cell. It uses the 
        center of mass (COM) or  the geometric center (GC) of each molecule 
        to determine the necessary translation.

        """
        # Retrieve the selected object
        vobject_id   = self.coordinates_combobox.get_vobject_id()
        
        system_id    = self.coordinates_combobox.get_system_id ()
        
        bb  = self.coordinates_combobox.get_active()
        
        print(' vobject_id', vobject_id, 
               '\n system_id',  system_id,
               '\n get_active', bb)
        #print(self.main.vm_session.vm_objects_dic.keys())
        
        vismol_object = self.main.vm_session.vm_objects_dic[vobject_id]
        
        
        # Get current selection (not directly used here, but may be relevant for extensions)
        selections = self.vm_session.selections[self.vm_session.current_selection]
        
        # Extract cell parameters (unit cell dimensions)
        a = vismol_object.cell_parameters['a']
        b = vismol_object.cell_parameters['b']
        c = vismol_object.cell_parameters['c']
        
        
        #---------------------------------------------------------------------------
        #                      Center the selection in the box
        #---------------------------------------------------------------------------
        _center_sel =  self.builder.get_object('chkbox_center_selection').get_active() 
        # move selection to the center of the box? 
        """
        In this step, we will move the geometric center of the selection 
        to the center of the box – all atoms are displaced, but the 
        selection serves as the reference.
        """
        if _center_sel:
            center_a = a/2
            center_b = b/2
            center_c = c/2
            
            
            for frame_index, frame in enumerate(vismol_object.frames):
                 
                x = 0.0 
                y = 0.0 
                z = 0.0 
                
                n = 0
                for atom in selections.selected_atoms:
                    xyz = atom.coords( frame = frame_index)
                    x += xyz[0]
                    y += xyz[1]
                    z += xyz[2]
                    n+=1
                
                x  = x/n # center of mass
                y  = y/n # center of mass
                z  = z/n # center of mass
                
                dx = x - center_a
                dy = y - center_b
                dz = z - center_c
                
                for atm_coords in frame:
                    atm_coords[0] += -dx
                    atm_coords[1] += -dy
                    atm_coords[2] += -dz
                    
                
        
        #'''
        #---------------------------------------------------------------------------
        #                  Put Everything Back Inside the Box
        #---------------------------------------------------------------------------
        # Iterate over all frames
        for frame_index, frame in enumerate(vismol_object.frames):
            
            #mol_idx = 0
            for mol_idx, molecule in enumerate(vismol_object.molecules.values()):
            #for molecule in vismol_object.molecules.values():
                x = 0.0 
                y = 0.0 
                z = 0.0 
                
                n = 0 
                
                indexes = []
                for atom in molecule.atoms.values():
                    xyz = atom.coords( frame = frame_index)
                    x += xyz[0]
                    y += xyz[1]
                    z += xyz[2]
                    indexes.append(atom.index -1)                
                    n+=1
                
                x = x/n #
                y = y/n # Geometric center
                z = z/n #
                
                #-----------------------------------------------------------
                '''
                X and O are outside the box;
                           ---------
                           |       |       O
                           |       |
                           |       |
                      X    ---------
               
                X is located before the origin of the space that defines
                the box. Therefore, in this case, we need to subtract 
                one (1) from the parameters na, nb, or nc (as applicable) 
                in order to correctly shift the molecules inside the 
                box in the next step.
                
                '''
                if x > 0:
                    na = int(x/a) # if x > a, so na = 0
                else:
                    na = int(x/a) - 1
                
                if y > 0:
                    nb = int(y/b)
                else:
                    nb = int(y/b)- 1
                
                if z > 0:
                    nc = int(z/c)
                else:
                    nc = int(z/c)- 1
                #-----------------------------------------------------------
        
                for index in indexes:
                    if x > a or x < 0 or y > b or y < 0 or z > c or z < 0:
                        vismol_object.frames[frame_index][index][0] += -(na*a)
                        vismol_object.frames[frame_index][index][1] += -(nb*b)
                        vismol_object.frames[frame_index][index][2] += -(nc*c)
        
        self.vm_session.set_frame(frame=0)
        self.vm_session.vm_glcore.queue_draw()    
        #'''

        
        



def get_subset (vismol_object, _type = 0, ignore_H = True, selections = None):
    """ Function doc """
    subset = []
   
   #.Protein
    if _type == 0:
        for atom in vismol_object.atoms.values():
            if atom.residue.is_protein:
                if atom.symbol == 'H':
                    pass
                else:
                    subset.append(atom.index -1)


    if _type == 1:
        for atom in vismol_object.atoms.values():
            if atom.residue.is_protein:
                if atom.name in  ["CA", 'C', 'O', 'N']:
                    subset.append(atom.index -1)
                    print(atom.name)

    
    if _type == 2:
        for atom in vismol_object.atoms.values():
            if atom.residue.is_protein:
                if atom.name in  ["CA", 'C', 'N']:
                    subset.append(atom.index -1)
                    print(atom.name)



    if _type == 3:
        for atom in vismol_object.atoms.values():
            if atom.residue.is_protein:
                if atom.name == "CA":
                    subset.append(atom.index -1)
       
    if _type == 4:
        if selections:
            for atom in selections.selected_atoms:
                if atom.symbol == 'H':
                    pass
                else:
                    subset.append(atom.index -1)
        pass
    return subset
                    
def align_trajectory(traj, reference=None, subset=None):
    """

    Aligns a trajectory to a reference frame using the Kabsch algorithm.
    -traj: np.ndarray with shape (n_frames, n_atoms, 3)
    -reference: frame used as reference (default = first frame)
    -subset: atom indices used for the alignment (e.g., backbone)
    
    
    
    kabsch_fixed(A, B)
        Computes the rotation R and translation t that minimize the RMSD 
        between points in B and A (with 1:1 correspondence).
        
        Ensures that the rotation matrix is a proper rotation 
        (determinant +1), correcting the reflection case after the SVD
        (Singular Value Decomposition).

    superpose_by_subset(A, B, subset_A, subset_B)
        
        Allows aligning the entire molecule B using only a subset of 
        atoms to compute R and t.
        
        This is typically used, for example, in molecular dynamics when
        one wants to align by the backbone or by the Cα atoms, but move 
        the whole structure.
    """
    
    if reference is None:
        reference = traj[0]

    aligned_traj = []
    for frame in traj:
        frame_aligned, R, t = superpose_by_subset(traj[reference], frame, subset, subset)
        aligned_traj.append(frame_aligned)

    return np.array(aligned_traj)

def kabsch_fixed(A, B):
    """
    A, B: (N,3) with a one-to-one correspondence (same N, same order).
    Returns R and t to map B → A.
    
    Parameters
    ----------
    A : array-like of shape (N, 3)
        Target coordinates. Points must correspond 1:1 with B and be in the same order.
    B : array-like of shape (N, 3)
        Source coordinates to be rotated/translated onto A.

    Returns
    -------
    R : ndarray of shape (3, 3)
        Proper rotation matrix (det(R) = +1).
    t : ndarray of shape (3,)
        Translation vector such that:  A ≈ (R @ B.T).T + t
    
    """
    # 1) Centroids of each point set (mean along the N points for x,y,z).
    cA = A.mean(axis=0)   # shape (3,)
    cB = B.mean(axis=0)   # shape (3,)

    # 2) Center both sets at the origin to remove translation before fitting rotation.
    A0 = A - cA           # shape (N, 3)
    B0 = B - cB           # shape (N, 3)

    # 3) Cross-covariance (or correlation) matrix between B and A.
    #    This captures how directions in B align with directions in A.
    H = B0.T @ A0         # shape (3, 3)

    # 4) Singular Value Decomposition of H: H = U * diag(S) * Vt.
    #    SVD provides orthonormal bases U and V that best align the two sets.
    U, S, Vt = np.linalg.svd(H)

    # 5) Optimal rotation that minimizes RMSD (Kabsch solution).
    #    Note: Vt.T is V, so R = V * U^T.
    R = Vt.T @ U.T        # shape (3, 3)

    # 6) Enforce a proper rotation (no reflection).
    #    If det(R) < 0, flip the sign of the last singular vector and recompute R.
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # 7) Translation that maps the centroid of B onto the centroid of A after rotation.
    #    Ensures B is positioned correctly relative to A.
    t = cA - R @ cB       # shape (3,)

    # 8) Return the rigid transform B → A.
    return R, t

def superpose_by_subset(A, B, subset_A=None, subset_B=None):
    """
    Compute a rigid transformation (rotation R and translation t) to align
    molecule B onto molecule A using only a subset of atoms, and then
    apply the transformation to the entire molecule B.

    Parameters
    ----------
    A : array-like of shape (N, 3)
        Target coordinates. Usually the reference frame.
    B : array-like of shape (N, 3)
        Source coordinates to be aligned onto A.
    subset_A : array-like of int, optional
        Indices of atoms in A to use for calculating the alignment.
        If None, all atoms are used.
    subset_B : array-like of int, optional
        Indices of atoms in B corresponding to subset_A.
        If None, all atoms are used.

    Returns
    -------
    B_aligned : ndarray of shape (N, 3)
        Coordinates of B after applying the rigid transformation.
    R : ndarray of shape (3, 3)
        Rotation matrix used for alignment.
    t : ndarray of shape (3,)
        Translation vector used for alignment.
    """

    # 1) If no subset is provided, use all atoms.
    if subset_A is None:
        subset_A = np.arange(A.shape[0])
    if subset_B is None:
        subset_B = np.arange(B.shape[0])

    # 2) Compute optimal rotation R and translation t using Kabsch
    #    on the selected subsets of atoms. This ensures that only
    #    the chosen “core” (e.g., backbone) is used for alignment.
    R, t = kabsch_fixed(A[subset_A], B[subset_B])

    # 3) Apply the computed transformation to the entire molecule B.
    #     @ -> In summary: @ is a cleaner and more modern way to perform 
    #     matrix multiplication in Python, without explicitly calling np.dot.
    B_aligned = (R @ B.T).T + t

    # 4) Return the aligned molecule and the transformation parameters.
    return B_aligned, R, t


