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
system = ImportSystem ('../pkl/bALA_OPLS.pkl')
# summary
system.Summary()

'''
Umbrella sampling using a simple distance reaction coordinate
 

'''

parameters = {'NmaxThreads': 2,
             
             'MD_parm': {'collisionFrequency'       : 25.0,
                         'folder'                   : '/home/fernando',
                         'integrator'               : 'Verlet',
                         'logFrequency'             : 10,
                         'normalDeviateGenerator'   : None,
                         'pressure'                 : 1.0,
                         'pressureControl'          : False,
                         'pressureCoupling'         : 2000.0,
                         'seed'                     : 44561,
                         'simulation_type'          : 'Umbrella_Sampling',
                         'steps_dc'                 : 5000,  # data collection
                         'steps_eq'                 : 2000,  # equilibrium
                         'temperature'              : 300.0,
                         'temperatureControl'       : True,
                         'temperatureCoupling'      : 0.1,
                         'temperatureScaleFrequency': 100,
                         'temperatureStart'         : 300.0,
                         'timeStep'                 : 0.001,
                         'trajectory_frequency'     : 10},
             
             
             
             
             'OPT_parm': {'maxIterations': 600,
                         'optimizer'     : 'ConjugatedGradient',
                         'rmsGradient'   : 0.1},
             
             
             'RC1': {'ATOMS'         : [5, 17],
                     'ATOM_NAMES'    : ['O5', 'H17'],
                     'dincre'        : 0.0,
                     'dminimum'      : 2.0595834656366816,
                     'force_constant': 20.0,
                     'nsteps'        : 10,
                     'rc_type'       : 'simple_distance',
                     'sigma_pk1pk3'  : None,
                     'sigma_pk3pk1'  : None},
             
             'RC2': None,
             
             'folder'           : '../scratch',
             'input_type'       : 1, # means = parallel
             'simulation_type'  : 'Umbrella_Sampling',
             'source_folder'    : '/home/fernando/programs/EasyHybrid3/examples/scripts/bALA_OPLS_PES_scan.ptGeo',
             'system'           : system,
             'traj_folder_name' : '1-BALA_OPLS_umb_sam',
             'vobj'             : None}



simObj =  pMethods.UmbrellaSampling()
simObj.run(parameters)



