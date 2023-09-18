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
This is the simplest definition of the reaction coordinate, based on the 
distance between two selected atoms, #1 and #2 

r = d(#1 - #2)

If the step size used is a positive value, the distance between atoms 
will increase at each scan step. In this case, the harmonic potential 
restraint applies only to the simple distance between atoms #1 and #2.

 
 
for more:

https://sites.google.com/view/easyhybrid/user-guide/pes-scans?authuser=1
 

'''

parameters = {'NmaxThreads'          : 1,                  # Use 1 for one-dimensional reaction coordinate scans 
                                                           
             'RC1': {'ATOMS'         : [0, 3],             # atom indices. For single-distance two atoms are needed
                     'ATOM_NAMES'    : ['O0', 'N3'],       # atom names.
                     
                     'dincre'        : 0.1,                # distance increment, use negative values to decrease the distance
                     'dminimum'      : 4.208733002805674,  # starting distance
                     'force_constant': 4000.0,             # Force constant used in harmonic restraint
                     'nsteps'        : 5,                  #
                     'rc_type'       : 'simple_distance',  # simple_distance / multiple-distance
                     
                     'sigma_pk1pk3'  : None,               # Required only when multiple-distance is used
                     'sigma_pk3pk1'  : None},              # Required only when multiple-distance is used
             
             'RC2': None,
             
             'folder'             : '../scratch',
             'initial_coordinates': None,
             
             'log_frequency'      : 50,
             'maxIterations'      : 600.0,                     #This is the maximum number of interactions allowed at each geometry optimization.
             'optimizer'          : 'ConjugatedGradient',
             'rmsGradient'        : 0.1,
             'simulation_type'    : None,
             'system'             : system,
             'traj_folder_name'   : '1-MOL_AM1_QC10_PES_scan_simple_distance', # rename as you wish / Trj out is a folder containg pkl files
             'traj_type'          : 'pklfolder',
             'vobject_name'       : None
             }

simObj =  pMethods.RelaxedSurfaceScan()
simObj.run(parameters)



