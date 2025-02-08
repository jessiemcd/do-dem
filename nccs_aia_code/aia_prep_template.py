from sys import path as sys_path
sys_path.append('/home/jmdunca2/jessie_notebooks/')

import aia_prep as ap

#place to save prepped map directories (aia_path)
map_save_path='/explore/nobackup/people/jmdunca2/'
#place to save DEM inputs
NCCS_save_path='./saved_AIA_dem_inputs/'
#Location of AIA temperature response
NCCS_aia_resp_path='/home/jmdunca2/aia_tresp_en.dat'
#Location of aia error tables
errortab='/home/jmdunca2/aia_V3_error_table.txt'

data_dir = 
f =


ap.file_prep(f, data_dir, NCCS_save_path, map_save_path,
                NCCS_aia_resp_path, errortab,
                clobber=False, aia_clobber=True)
