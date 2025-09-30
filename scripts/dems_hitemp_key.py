#Path to top-level do-dem directory - edit for your system.
path_to_dodem = '/Users/jmdunca2/do-dem/'
from sys import path as sys_path
sys_path.append(path_to_dodem+'/dodem/')

import all_nu_analysis as ana
import pickle
import importlib

key = 


#rerun all, (hightemp analysis to do rescale)

importlib.reload(ana)

with open('all_targets.pickle', 'rb') as f:
	data = pickle.load(f)

if key =='22-apr-16_1':
	print('22-apr-16_1 failed TIS, nothing to do.')
else:
	ARDict = data[key]
	id_dirs = ARDict['datapaths']
	working_dir = ARDict['working_dir']
	method = ARDict['method']
	
	ARDict['prepped_aia'] = working_dir+'all_aia_dicts_'+key+'_post/'
	
	data[key] = ARDict
	
	with open('all_targets.pickle', 'wb') as f:
		 # Pickle the 'data' dictionary using the highest protocol available.
		 pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 
	
	
	ana.do_key_dem(key, plot_xrt=False, method=method, high_temp_analysis=True, rscl=True) #missing_last=True, missing_orbit=1, plot_xrt=False)
	print('')
	print('')

print('Finished!')