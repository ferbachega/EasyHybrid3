import sys
from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                  
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                 
from pSimulation               import *

sys.path.insert(1, '/home/fernando/programs/EasyHybrid3/pdynamo')
import p_methods  as pMethods



# importing system
system = ImportSystem ('../pkl/molecules_AM1.pkl')
# summary
system.Summary()



'''
This reaction coordinate is defined by a combination of distances 
obtained from the selection of three atoms, #1, #2, and #3, as shown 
in Figure 3. In this case, the reaction coordinate is:

r = d(#1 - #2) - d(#2 - #3)

where d(#1 - #2) and d(#2 - #3) are the distances between atoms #1 and 
#2, and #2 and #3, respectively. If the step size is a positive value, 
the distance between #1 and #2 decreases, and the distance between atoms 
#2 and #3 increases, resulting in the tendency for the bond between #1 
and #2 to break, and the formation of the bond between #2 and #3. In 
this case, the harmonic potential is restrictive for both distances that 
constitute the reaction coordinate simultaneously.

It is important to emphasize that this does not characterize a 
two-dimensional scan, where each distance comprises a dimension, but 
rather, when dragging the atoms in the defined reaction coordinate, they
 will be restricted by the two defined distances simultaneously.
 
 
for more:

https://sites.google.com/view/easyhybrid/user-guide/pes-scans?authuser=1
 
'''

parameters = {'NmaxThreads'          : 1,                   # Use 1 for one-dimensional reaction coordinate scans 
                                                           
             'RC1': {'ATOMS'         : [7, 9, 3],           # atom indices. For single-distance two atoms are needed
                     'ATOM_NAMES'    : ['O7', 'H9', 'N3'],  # atom names.
                     
                     'dincre'        : 0.1,                # distance increment
                     'dminimum'      : -0.713144174206749, # starting distance
                     'force_constant': 4000.0,             # Force constant used in harmonic restraint
                     'nsteps'        : 10,                 #
                     'rc_type'       : 'multiple_distance',  # simple_distance / multiple-distance
                     
                     'sigma_pk1pk3'  :  0.5332049150006165,   # Required only when multiple-distance is used
                     'sigma_pk3pk1'  : -0.466795084999383},   # Required only when multiple-distance is used
             
             'RC2': None,
             
             'folder'             : '../scratch',
             'initial_coordinates': None,
             
             'log_frequency'      : 50,
             'maxIterations'      : 600.0,
             'optimizer'          : 'ConjugatedGradient',
             'rmsGradient'        : 0.1,
             'simulation_type'    : None,
             'system'             : system,
             'traj_folder_name'   : '1-MOL_AM1_QC10_PES_scan_multiple_distance', # rename as you wish / Trj out is a folder containg pkl files
             'traj_type'          : 'pklfolder',
             'vobject_name'       : None
             }

simObj =  pMethods.RelaxedSurfaceScan()
simObj.run(parameters)



