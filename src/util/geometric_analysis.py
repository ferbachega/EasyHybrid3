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
import numpy as np


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
    

def get_dihedral(vobject, index1, index2, index3, index4):
    # Convert the coordinates to numpy arrays
    #atom1 = vobject.atoms[index1]
    #atom2 = vobject.atoms[index2]
    #atom3 = vobject.atoms[index3]
    #atom4 = vobject.atoms[index4]
    
    a1_coord = vobject.atoms[index1].coords()
    a2_coord = vobject.atoms[index2].coords()
    a3_coord = vobject.atoms[index3].coords()
    a4_coord = vobject.atoms[index4].coords()
    
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




