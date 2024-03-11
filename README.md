# do-dem
Code to prepare data from NuSTAR and other instruments for use in DEM analysis, do DEMs (either using DEMREG or xrt_dem_iterative2), and visualize results.

For specifics on what was done in the May 29 2018 AR DEM analysis, see the Duncan_et_al_2024_Process folder.

For information on how to get this going on your own machine, see below.

Requirements:
- Python packages: numpy, matplotlib, sunpy, astropy, aiapy, scipy, regions
- NuSTAR pysolar: https://github.com/ianan/nustar_pysolar/tree/main 
- For NuSTAR data prep: Heasoft with NuSTAR mission-specific tools, installed and initialized (https://heasarc.gsfc.nasa.gov/docs/software/heasoft/). 
- SSWIDL installation (https://www.lmsal.com/solarsoft/ssw_setup.html) with all relevant observatories included. This code explicitly uses the SDO/AIA, Hinode/XRT, HESSI, and XRAY libraries (see the hissw section below). In addition, SSW provides:
  - AIA error tables (set path to error tables in your system in aia_dem_prep.py)
  - Functionality for XRT data download and prep to level 1.
- hissw (https://wtbarnes.github.io/hissw/) – Package used to run SSWIDL from within Python. Make sure to set up your hisswrc file (see installation instructions). This is used for:
  - fetching aia response (using aia_get_response.pro) - done in aia_dem_prep.py
  - fetching xrt response (using make_xrt_temp_resp.pro, make_xrt_wave_resp.pro, etc) - done in xrt_dem_prep.py
  - fetching x-ray thermal emission model for use in generating NuSTAR temperature response (using fvth.pro) - done in nustar_dem_prep.py
  - doing DEMs via xrt_dem_iterative2.pro (wrapper for this method included in dodem.py)
- To use the demreg method for DEM calculation (python version): demreg repository (https://github.com/ianan/demreg/tree/master) Note that this has its own requirements, which are described in https://github.com/ianan/demreg/tree/master/python 


Procedure:

- Download NuSTAR data from HEASARC https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3browse.pl
  - Unpack files in data directories (bash: gunzip -r (path to OBSID-named directory))
  - Run initial pipeline (with pixpos correction) - if you don't know how to do this, contact Jessie or any other NuSTAR Heliophysics Team member. 
- If desired, search for, download, and prep XRT data to level 1.
  - IDL procedure lets_get_this_xrt.pro (in other_idl) will allow you to do this via a search for XRT images at the same time and with an interesecting FOV with NuSTAR. Make sure to check the documentation for information about needed databases through SSWDB.
  - Or, do it your own way – information can be found in https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf. You will need level 1 data, and the XRT grade maps as well as the data files.
 
Examples:
- **getting_started.ipynb**: Initial DEM process example notebook (examples with and without XRT, with and without NuSTAR pile-up correction)
- **getting_started-xrt_iterative.ipynb**: Similar, but using xrt_dem_iterative2.pro instead of DEMREG
- **time_interval_selection.ipynb**: Example of how to select time intervals within a longer interval (or entire NuSTAR orbit) in order to do time-resolved DEM analysis while ensuring sufficient NuSTAR statistics in your highest energy range. 
- COMING SOON: examples of how to do... variable temperature range DEMs, DEMs where you insert your own PHA files (like if you modify them to remove a non-thermal component...), others - have a request? Contact Jessie.
