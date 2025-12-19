#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
#
#  Copyright 2022-2025 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
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
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      Provides functions for selecting atoms and residues in pDynamo systems
#      to facilitate QM/MM partitioning and molecular simulations.
#
from collections import deque
import numpy as np
import math

def get_simple_distance (a1_coord, a2_coord):
    """ Function doc """ 
    dx = a1_coord[0] - a2_coord[0]
    dy = a1_coord[1] - a2_coord[1]
    dz = a1_coord[2] - a2_coord[2]
    dist = (dx**2+dy**2+dz**2)**0.5
    #print('distance a1 - a2:', dist)
    return dist


def get_distance (vobject, index1, index2):
    """ Function doc """
    #print( index1, index2)
    atom1 = vobject.atoms[index1]
    atom2 = vobject.atoms[index2]
    a1_coord = atom1.coords()
    a2_coord = atom2.coords()
    
    dx = a1_coord[0] - a2_coord[0]
    dy = a1_coord[1] - a2_coord[1]
    dz = a1_coord[2] - a2_coord[2]
    dist = (dx**2+dy**2+dz**2)**0.5
    #print('distance a1 - a2:', dist)
    return dist


def get_simple_angle(a1_coord, a2_coord, a3_coord):
    # Convert the coordinates to numpy arrays
    atom1 = np.array(a1_coord)
    atom2 = np.array(a2_coord)
    atom3 = np.array(a3_coord)
    
    # Compute vectors between the atoms
    vec1 = atom1 - atom2
    vec2 = atom3 - atom2
    
    # Normalize the vectors
    vec1 /= np.linalg.norm(vec1)
    vec2 /= np.linalg.norm(vec2)
    
    # Compute the dot product
    dot = np.dot(vec1, vec2)
    
    # Compute the angle
    angle = np.arccos(dot)
    angle = np.degrees(angle)
    
    return angle 
    

def get_angle (vobject, index1, index2, index3):
    
    atom1 = vobject.atoms[index1]
    atom2 = vobject.atoms[index2]
    atom3 = vobject.atoms[index3]

    a1_coord = atom1.coords()
    a2_coord = atom2.coords()
    a3_coord = atom1.coords()
 
    # Convert the coordinates to numpy arrays
    atom1 = np.array(a1_coord)
    atom2 = np.array(a2_coord)
    atom3 = np.array(a3_coord)

    # Compute vectors between the atoms
    vec1 = atom1 - atom2
    vec2 = atom3 - atom2

    # Normalize the vectors
    vec1 /= np.linalg.norm(vec1)
    vec2 /= np.linalg.norm(vec2)

    # Compute the dot product
    dot = np.dot(vec1, vec2)

    # Compute the angle
    angle = np.arccos(dot)
    angle = np.degrees(angle)

    return angle 


def get_simple_dihedral(a1_coord, a2_coord, a3_coord, a4_coord):
    """ Function doc """
    atom1 = np.array(a1_coord)
    atom2 = np.array(a2_coord)
    atom3 = np.array(a3_coord)
    atom4 = np.array(a4_coord)
    
    #print(atom1, atom2,atom3  , atom4)
    
    # Compute vectors between the atoms
    vec1 = atom2 - atom1
    vec2 = atom3 - atom2
    vec3 = atom4 - atom3
    
    # Normalize the vectors
    vec1 /= np.linalg.norm(vec1)
    vec2 /= np.linalg.norm(vec2)
    vec3 /= np.linalg.norm(vec3)
    
    # Compute the cross products
    cross1 = np.cross(vec1, vec2)
    cross2 = np.cross(vec2, vec3)
    
    # Compute the dot product between the cross products
    dot = np.dot(cross1, cross2)
    
    # Compute the dihedral angle
    dihedral = np.arctan2(np.linalg.norm(np.cross(cross1, cross2)), dot)
    dihedral = np.degrees(dihedral)
    return dihedral
    

#def get_dihedral(vobject, index1, index2, index3, index4):
def get_dihedral(atom1, atom2, atom3, atom4):
    # Convert the coordinates to numpy arrays
    #atom1 = vobject.atoms[index1]
    #atom2 = vobject.atoms[index2]
    #atom3 = vobject.atoms[index3]
    #atom4 = vobject.atoms[index4]
    
    a1_coord = atom1.coords()
    a2_coord = atom2.coords()
    a3_coord = atom3.coords()
    a4_coord = atom4.coords()
    
    #print(a1_coord, a2_coord,a3_coord  , a4_coord)
    
    atom1 = np.array(a1_coord)
    atom2 = np.array(a2_coord)
    atom3 = np.array(a3_coord)
    atom4 = np.array(a4_coord)

    #print(atom1, atom2,atom3  , atom4)

    # Compute vectors between the atoms
    vec1 = atom2 - atom1
    vec2 = atom3 - atom2
    vec3 = atom4 - atom3

    # Normalize the vectors
    vec1 /= np.linalg.norm(vec1)
    vec2 /= np.linalg.norm(vec2)
    vec3 /= np.linalg.norm(vec3)

    # Compute the cross products
    cross1 = np.cross(vec1, vec2)
    cross2 = np.cross(vec2, vec3)

    # Compute the dot product between the cross products
    dot = np.dot(cross1, cross2)

    # Compute the dihedral angle
    dihedral = np.arctan2(np.linalg.norm(np.cross(cross1, cross2)), dot)
    dihedral = np.degrees(dihedral)
    return dihedral







def center_on_atom(vobject, x, y, z, subgroup, frame = 0):
    '''the function that centers on the x, y and z coordinates  '''
    #print(x, y, z)
    #for index in range(len(vobject.atoms.keys()) ):
    for index in subgroup:
        vobject.frames[frame][index][0] = vobject.frames[frame][index][0] - x
        vobject.frames[frame][index][1] = vobject.frames[frame][index][1] - y
        vobject.frames[frame][index][2] = vobject.frames[frame][index][2] - z
    return vobject

def move_to_origin (molecule, center_id):
    """ Function doc """
    coord_center = molecule[center_id] 
    print(coord_center)
    #coord_center = ['x',10,10,10] 
    
    for index, atom in enumerate(molecule):
        print(index, atom)
        molecule[index][1] = molecule[index][1] - coord_center[1] 
        molecule[index][2] = molecule[index][2] - coord_center[2] 
        molecule[index][3] = molecule[index][3] - coord_center[3] 
    
    return molecule

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation
    about
    the given axis by theta radians.
    """
    axis  = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2)
    b, c, d = -axis * math.sin(theta / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def rotate_dihedral(vismol_object = None       , 
                         subgroup = None       , 
                            #center= ['',0,0,0], 
                             axis = None       , 
                            theta = 0.500      ,
                            frame = 0 ):
    """ Function doc """
    for index  in subgroup:
        pos  = (vismol_object.frames[frame][index][0],
                vismol_object.frames[frame][index][1],
                vismol_object.frames[frame][index][2],
                )
           
        
        new_pos = np.dot(rotation_matrix(axis,
                                         theta),
                                         pos
                                         )
                                         
        vismol_object.frames[frame][index][0] = new_pos[0]
        vismol_object.frames[frame][index][1] = new_pos[1]
        vismol_object.frames[frame][index][2] = new_pos[2]
    return vismol_object

def rotate_bond (vobject     =  None, 
                index_center =  None, 
                index_vector =  None, 
                subgroup     =  None, 
                theta        =  1, 
                frame        =  0 ):
    """ Function doc """
    
    #print (vobject.frames[frame][index_center])
    #print (vobject.frames[frame][index_vector])
    
    
    #print (vobject.frames[0][index_center])
    #print (vobject.frames[0][index_center])
    
    center = (vobject.frames[frame][index_center][0],
              vobject.frames[frame][index_center][1],
              vobject.frames[frame][index_center][2])
    #      
    center_on_atom(vobject,
                   center[0],
                   center[1],
                   center[2],
                   subgroup,
                   frame)
    
    axis = (
            vobject.frames[frame][index_vector][0],                            
            vobject.frames[frame][index_vector][1],                            
            vobject.frames[frame][index_vector][2]                             
            )
    
    vobject = rotate_dihedral(vismol_object   = vobject, 
                              subgroup = subgroup,
                              #center   = center,  
                              axis = axis, 
                              theta = theta,
                              frame = frame)
    
    vobject = center_on_atom(vobject,  center[0]*-1,
                                       center[1]*-1,
                                       center[2]*-1,
                                       subgroup,
                                       frame)    

def rotate_bond_from_indexes (molecule, index_center, index_vector, theta):
    """ Function doc """
    subgroup = find_subgroup (index_center, index_vector)
    molecule = rotate_bond (molecule, index_center, index_vector, subgroup, theta )
    return molecule

def find_subgroup(atom1, atom2, top):
    
    #Versão iterativa, eficiente e segura contra ciclos.
    #Retorna todos os átomos acessíveis a partir de atom2,
    #excluindo atom1 e evitando loops cíclicos.
    #
    #Parâmetros:
    #    atom1 : int
    #        Átomo "pai" inicial — usado para evitar travessia de volta.
    #
    #    atom2 : int
    #        Átomo inicial da busca (onde deve-se começar).
    #
    #    top : dict[int, list[int]]
    #        Conectividade do sistema (grafo). Exemplo:
    #        top[i] = [lista de átomos conectados ao átomo i]
    #
    #Retorno:
    #    list[int]  — lista dos átomos pertencentes ao subgrupo.
    """
    Iterative, efficient, and cycle-safe version.
    Returns all atoms accessible from atom2, excluding atom1 and avoiding cyclic loops.

    Parameters:

    atom1 : int
    Initial “parent” atom — used to prevent traversal back.

    atom2 : int
    Starting atom for the search (the traversal origin).

    top : dict[int, list[int]]
    System connectivity (graph). Example:
    top[i] = [list of atoms connected to atom i]

    Return:

    list[int] — list of atoms belonging to the subgroup.
    """

    subgroup = []
    visited = set([atom1])   # marca atom1 como visitado para nunca voltar
    queue = deque([atom2])   # busca em profundidade ou largura (DF/ BF)

    while queue:
        current = queue.pop()
        
        for neighbor in top[current]:
            if neighbor in visited:
                continue
            
            visited.add(neighbor)
            subgroup.append(neighbor)

            # só avança se o átomo tiver ramificações
            if len(top[neighbor]) > 1:
                queue.append(neighbor)

    return subgroup




