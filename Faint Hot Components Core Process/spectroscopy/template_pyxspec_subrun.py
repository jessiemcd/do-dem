import subprocess
import os

pathstring = '/Users/jmdunca2/do-dem/initial_dem_10oct17/region_0/01-38-30_01-40-15/'


isoscript = pathstring+'/isot_pyxspec.py'
isoout = pathstring+'/isot_out.txt'

twotscript = pathstring+'/twot_pyxspec.py'
twotout = pathstring+'/twot_out.txt'

tntscript = pathstring+'/tnt_pyxspec.py'
tntout = pathstring+'/tnt_out.txt'

try:
   result = subprocess.run('python '+isoscript+' > '+isoout, shell=True, capture_output=False, timeout=600)
except subprocess.TimeoutExpired:
   print('Timed out on Isothermal.')
   print('Path: ', pathstring)

try:
   result = subprocess.run('python '+twotscript+' > '+twotout, shell=True, capture_output=False, timeout=600)
except subprocess.TimeoutExpired:
   print('Timed out on Double Thermal.')
   print('Path: ', pathstring)

try:
    result = subprocess.run('python '+tntscript+' > '+tntout, shell=True, capture_output=False, timeout=600)
except subprocess.TimeoutExpired:
    print('Timed out on Thermal + Broken Power Law.')
    print('Path: ', pathstring)