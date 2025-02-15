#Path to top-level do-dem directory - edit for your system.
path_to_dodem = '/Users/jmdunca2/do-dem/'
from sys import path as sys_path
sys_path.append(path_to_dodem+'/dodem/')

import time_interval_selection as tis

import pickle

key = "20-jul-21"

with open('all_targets.pickle', 'rb') as f:
    data = pickle.load(f)

index = 3
tis.one_orbit_tis_wrapper(key, data, index, method='doublegauss', use_set_regionfiles=True)
