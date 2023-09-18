#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  install.py
#  
#  Copyright 2021 Fernando <fernando@winter>
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
import os



welcome_text = '''
Welcome to the EasyHydrid3 installation! We're thrilled that you've chosen to 
install our software and we hope that it will be a valuable addition to your 
propose. 

EasyHybrid is an open-source graphical environment that is specifically 
designed for the pDynamo3 package. The software is developed using Python3, 
Cython3, GTK, and modern OpenGL, and incorporates the VisMol graphical engine 
to render 3D structures. As an extension of pDynamo, EasyHybrid provides a 
user-friendly interface that allows users to execute fundamental pDynamo 
routines. It also provides an efficient way to inspect, edit, and export 
pDynamo systems for further simulations by leveraging Python scripting in text 
mode. EasyHybrid3 has been released under the GPL3 license.
'''



final_text = '''
 - - - Congratulations! The installation of EasyHydrid3 was successful. - - -  

Don't forget to source the environment_bash.com file before start. You can also
add the environment_bash.com path in your .bashrc inicialization file, located 
in your home directory.. 


Thank you for choosing EasyHydrid and Happy Simulation!
'''

#Once again, thank you for choosing [Software Name] and we hope you enjoy using our 
#software!
#'''

print(welcome_text)


print('Checking for pDynamo3 installation:')

try:
    from pBabel                    import *                                     
    print('pBabel...............OK!')
except:
    print('pBabel...............FAILED!')

try:    
    from pCore                     import *
    print('pCore................OK!')
except:
    print('pCore................FAILED!')

try:    
    from pMolecule                 import *                  
    print('pMolecule............OK!')

except:
    print('pMolecule............FAILED!')

try:    
    from pScientific               import *                                     
    print('pScientific..........OK!')

except:
    print('pScientific..........FAILED!')


try:    
    from pSimulation               import *
    print('pSimulation..........OK!')
except:
    print('pSimulation..........FAILED!')







print('\nChecking for OpenGL installation:')
try:
    from OpenGL import GL
          #pSimulation..........FAILED!
    print('OpenGL...............OK!')
except:
    print('OpenGL...............FAILED!')




print('\nChecking for numpy installation:')
try:
    import numpy as np
    print('numpy................OK!')
except:
    print('numpy................FAILED!')




print('\nWriting the file: environment_bash.com')

try:
    #os.system('python3 setup.py build_ext --inplace')
    envfile = open('environment_bash.com', 'w')
    env     = os.getcwd()

    text = ''
    text += '#!/bin/bash\n'
    text += '''# . Bash environment variables and paths to be added to a user's ".bash_profile" file.\n'''
    text += '# . The root of the program.\n'
    text += 'EASYHYBRID_HOME={} ; export EASYHYBRID_HOME\n'.format(env)
    #text += 'PYTHONPATH=$VISMOL_HOME ; export PYTHONPATH\n'
    envfile.write(text)
    print (final_text)




except:
    print('Failed! File "environment_bash" cannot be written.')







