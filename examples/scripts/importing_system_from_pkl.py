#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  importing_system_from_pkl.py
#  
#  Copyright 2023 Fernando Bachega <fernando@fernando-Inspiron-7537>
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


'''
This is a basic script to open an already prepared system in pDynamo. 
It simply consists of loading a pkl file and doing a quick system 
assessment.
'''

from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                  
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                 
from pSimulation               import *


# importing system
system = ImportSystem ('../pkl/bALA_eq_OPLS.pkl')

# summary
system.Summary()

