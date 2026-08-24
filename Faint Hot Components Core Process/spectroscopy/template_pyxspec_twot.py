#Path to top-level do-dem directory - edit for your system.
path_to_dodem = '/Users/jmdunca2/do-dem/'
from sys import path as sys_path
sys_path.append(path_to_dodem+'/dodem/')

import pyxspec_bb as pxb
import pathlib

pathstring = '/Users/jmdunca2/do-dem/initial_dem_10oct17/region_0/01-38-30_01-40-15/'
plot=True


pha_path = pathlib.Path(pathstring)


archive, model_name = pxb.do_xspectral_fit(pha_path, expression='const*(vapec+vapec)',
                           model_name = 'two_thermal_with_pileup_nogain',
                           fit_slope=False,
                           plot=plot)

# gain_slope_a = getattr(archive.instruments['FPM A'], model_name).response_parameters['slope'].value
# gain_slope_b = getattr(archive.instruments['FPM B'], model_name).response_parameters['slope'].value
# 
# print(gain_slope_a, gain_slope_b)
# 
# print('')
# print('')
# 
# slopes = {'FPM A': gain_slope_a,
#           'FPM B': gain_slope_b
#          }
# 
# import pickle
# 
# with open(pathstring+'/xspec_out/gain_slopes.pickle', 'wb') as f:
#     # Pickle the 'data' dictionary using the highest protocol available.
#     pickle.dump(slopes, f, pickle.HIGHEST_PROTOCOL)