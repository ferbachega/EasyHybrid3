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
                    This is  2D SCAN!
              
              
              
              
                 RC1| 
                    |
                    |
                    |
                    |
                    |
                    |
                    |
                    |-----------------
                                    RC2


The first reaction coordinate is sampled by a simple distance constraint. 

The second reaction coordinate is sampled using a multiple distance 
constraint.

In this case most of the process can be done in parallel (NmaxThreads > 1)

for more:

https://sites.google.com/view/easyhybrid/user-guide/pes-scans?authuser=1
 

'''

parameters = {'NmaxThreads'          : 4,                  # Use 1 for one-dimensional reaction coordinate scans 
                                                           
             'RC1': {'ATOMS'         : [0, 3],             # atom indices. For single-distance two atoms are needed
                     'ATOM_NAMES'    : ['O0', 'N3'],       # atom names.
                     
                     'dincre'        : 0.1,                # distance increment, use negative values to decrease the distance
                     'dminimum'      : 4.208733002805674,  # starting distance
                     'force_constant': 4000.0,             # Force constant used in harmonic restraint
                     'nsteps'        : 5,                  #
                     'rc_type'       : 'simple_distance',  # simple_distance / multiple-distance
                     
                     'sigma_pk1pk3'  : None,               # Required only when multiple-distance is used
                     'sigma_pk3pk1'  : None},              # Required only when multiple-distance is used
             
             
             'RC2': {'ATOMS'         : [7, 9, 3],          
                     'ATOM_NAMES'    : ['O7', 'H9', 'N3'], 
                     
                     'dincre'        : 0.1,                
                     'dminimum'      : -0.713144174206749, 
                     'force_constant': 4000.0,             
                     'nsteps'        : 10,                 
                     'rc_type'       : 'multiple_distance',
                     
                     'sigma_pk1pk3'  :  0.5332049150006165,
                     'sigma_pk3pk1'  : -0.466795084999383},
             
             
             'folder'             : '../scratch',
             'initial_coordinates': None,
             
             'log_frequency'      : 50,
             'maxIterations'      : 600.0,                     #This is the maximum number of interactions allowed at each geometry optimization.
             'optimizer'          : 'ConjugatedGradient',
             'rmsGradient'        : 0.1,
             'simulation_type'    : None,
             'system'             : system,
             'traj_folder_name'   : '1-MOL_AM1_QC10_PES_scan2D', # rename as you wish / Trj out is a folder containg pkl files
             'traj_type'          : 'pklfolder',
             'vobject_name'       : None
             }

simObj =  pMethods.RelaxedSurfaceScan()
simObj.run(parameters)



