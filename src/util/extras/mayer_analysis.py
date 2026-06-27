#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  script.py
#  
#  Copyright 2026 fernando <fernando@bahamuth>
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
import re

files  = os.listdir('.')

data = {}


def read_data (_file):
    """ Function doc 
 ATOM       NA         ZA         QA         VA         BVA        FA
  0 C      6.3232     6.0000    -0.3232     3.8657     3.8657    -0.0000
  1 H      0.8419     1.0000     0.1581     0.9580     0.9580    -0.0000
  2 H      0.8611     1.0000     0.1389     0.9842     0.9842    -0.0000
  3 C      5.8721     6.0000     0.1279     3.6664     3.6664    -0.0000
    
    """
    rawdata = open(_file, 'r')
    lines = rawdata.readlines()
    #print (rawdata)
    
    inside_block = False
    mpop_block = False
    mbond_block = False
    lchrg_block = False
    
    fdata = {
             'symbol'  : [],
             'lcharges': [],
             'NA'      : [],
             'ZA'      : [],
             'QA'      : [],
             'VA'      : [],
             'BVA'     : [],
             'FA'      : [],
             'mbo'     : [],
             }
    
    for i, line in enumerate(lines):
        
        #LOEWDIN ATOMIC CHARGES
        if "LOEWDIN ATOMIC CHARGES" in line:
            lchrg_block = True
            lchrg_start = i
        
        if 'LOEWDIN REDUCED ORBITAL CHARGES' in line:
            #Lchrg_block = True
            lchrg_end = i        
        
        #print(line)
        if "MAYER POPULATION ANALYSIS" in line:
            mpop_block = True
            mpop_start = i
            continue
        
        
        if 'Mayer bond orders larger than' in line:
            mbond_block = True
            mbond_start = i
            
            mpop_end = i
        
        if 'Environment variable NBOEXE' in line:
            #mbond_block = True
            mbond_end = i
    
    
    if lchrg_block:
        for j in range(lchrg_start, lchrg_end):   
            #print (lines[j].strip())
            l = lines[j].split()
            
            if len(l ) == 4:
                chrg = l[-1]
                fdata['lcharges'].append(float(chrg))
                fdata['symbol'].append(l[1])

        
    #'''    
    if mpop_block:
        for j in range(mpop_start, mpop_end):   
            #print (lines[j].strip())
            l = lines[j].split()
            if len(l ) ==8:
                fdata['NA' ].append(float(l[2]))
                fdata['ZA' ].append(float(l[3]))
                fdata['QA' ].append(float(l[4]))
                fdata['VA' ].append(float(l[5]))
                fdata['BVA'].append(float(l[6]))
                fdata['FA' ].append(float(l[7]))
                
                
                
                #print (lines[j])   
  
    if mbond_block:
        for j in range(mbond_start+1, mbond_end):   
            #print (lines[j].strip())
            l = lines[j].split('B')
            #if len(l ) ==8:
            
            for bond in l:
                if bond =='':
                    pass
                else:
                    fdata['mbo'].append(bond)
            

    return fdata


folder = '/home/fernando/programs/EasyHybrid3_dev/workspace/CCL_BR/6-jidM_XTB QC MODEL_QC6_umb_sam/ref/window0'
os.chdir(folder)

for _file in files:
    if _file[-3:] == 'log':
        fsplit = _file.split('.')
        
        try:
            fnum = fsplit[0].replace('frame', '')
            fnum = int(fnum)
            data[fnum] = {'file': _file}
            
            #print(_file)
            fdata = read_data (_file)
            
            data[fnum] = {'data': fdata} 
            #print(fdata)
        except:
            pass


#prop = 'lcharges'
prop = 'VA'
prop = 'QA'

#keys = data.keys()
stringline = ''
lines = []
for n , symbol in enumerate(data[0]['data']['symbol']):
    stringline = ''
    
    stringline += str(n)+' '+symbol+' '
    for i in range(len(data)):
    
        stringline+= str(data[i]['data'][prop][n])+ ' '
    
    #stringline+='\n'
    lines.append(stringline)
    print(stringline)
#print(lines)
        #for j in data[i]['data']['NA'] :
        
    


#
#
#
#_max = max(data.keys())
#
#for i in range(_max):
#    pass



#print(data)

