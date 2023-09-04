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
system = ImportSystem ('../pkl/waterTrimer_OPLS.pkl')
# summary
system.Summary()

'''
Umbrella sampling using a simple distance reaction coordinate
 

'''

parameters = {'MD_parm': {'collisionFrequency': 25.0,
                         'folder': '/home/fernando',
                         'integrator': 'Verlet',
                         'logFrequency': 10,
                         'normalDeviateGenerator': None,
                         'pressure': 1.0,
                         'pressureControl': False,
                         'pressureCoupling': 2000.0,
                         'seed': 44561,
                         'simulation_type': 'Umbrella_Sampling',
                         'steps_dc': 5000,
                         'steps_eq': 2000,
                         'temperature': 300.0,
                         'temperatureControl': True,
                         'temperatureCoupling': 0.1,
                         'temperatureScaleFrequency': 100,
                         'temperatureStart': 300.0,
                         'timeStep': 0.001,
                         'trajectory_frequency': 10},
             'NmaxThreads': 4,
             
             'OPT_parm': {'maxIterations': 600,
                         'optimizer': 'ConjugatedGradient',
                         'rmsGradient': 0.1},
             
             'RC1': {'ATOMS': [0, 6],
                     'ATOM_NAMES': ['O0', 'O6'],
                     'dincre': 0.0,
                     'dminimum': 2.835460117056201,
                     'force_constant': 20.0,
                     'nsteps': 10,
                     'rc_type': 'simple_distance',
                     'sigma_pk1pk3': None,
                     'sigma_pk3pk1': None},
             
             'RC2': {'ATOMS': [3, 0],
                     'ATOM_NAMES': ['O3', 'O0'],
                     'dincre': 0.0,
                     'dminimum': 2.838303644175005,
                     'force_constant': 20.0,
                     'nsteps': 10,
                     'rc_type': 'simple_distance',
                     'sigma_pk1pk3': None,
                     'sigma_pk3pk1': None},
             
             'folder': '../scratch',
             'input_type': 1,
             'simulation_type': 'Umbrella_Sampling',
             'source_folder': '../scripts/0-MolSys_OPLS_PES_scan.ptGeo',
             'system': system,
             'traj_folder_name': '1-WT_OPLS_umb_sam',
             'vobj': None}



simObj =  pMethods.UmbrellaSampling()
simObj.run(parameters)



