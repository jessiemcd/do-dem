# dodem
Code to prepare data from NuSTAR and other instruments for use in DEM analysis, do DEMs (either using DEMREG or xrt_dem_iterative2), and visualize results.

Requirements:
- Python packages: numpy, matplotlib, sunpy, astropy, aiapy, scipy
- For NuSTAR data prep: Heasoft with NuSTAR mission-specific tools, installed and initialized (https://heasarc.gsfc.nasa.gov/docs/software/heasoft/). 
- SSWIDL installation (https://www.lmsal.com/solarsoft/ssw_setup.html) with all relevant observatories included. This code explicitly uses the SDO/AIA and Hinode/XRT libraries (see the hissw section below). In addition, SSW provides:
  - AIA error tables (set path to error tables in your system in aia_dem_prep.py)
  - Functionality for XRT data download and prep to level 1.
- hissw: https://wtbarnes.github.io/hissw/ Package used to run SSWIDL from within Python. This is used for:
  - fetching aia response (using aia_get_response.pro) - done in aia_dem_prep.py
  - fetching xrt response (using make_xrt_temp_resp.pro, make_xrt_wave_resp.pro, etc) - done in xrt_dem_prep.py
  - fetching x-ray thermal emission model for use in generating NuSTAR temperature response (using fvth.pro) - done in nustar_dem_prep.py
  - doing DEMs via xrt_dem_iterative2.pro (wrapper for this method included in dodem.py)

