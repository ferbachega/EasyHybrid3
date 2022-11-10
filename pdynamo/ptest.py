from pDynamo2Easyhybrid import *
ps = pDynamoSession()

ps.load_a_new_pDynamo_system_from_dict (input_files = {'amber_prmtop': '/home/fernando/programs/VisMol/examples/TIM/7tim.top', 
                                                       'coordinates': '/home/fernando/programs/VisMol/examples/TIM/7tim.crd'}, 
                                                        system_type = 0, 
                                                        name = 'UNK', tag = 'tag')
