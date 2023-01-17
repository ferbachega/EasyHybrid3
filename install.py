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


#ef main(args):
#   return 0
#
#f __name__ == '__main__':
#   import sys
#   sys.exit(main(sys.argv))
