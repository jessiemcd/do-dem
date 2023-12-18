# This needs nustar_pysolar to be installed https://github.com/NuSTAR/nustar_pysolar
from sys import path as sys_path

# Change to your local copy's location...
sys_path.append('/Users/jessieduncan/demreg/python/')
from dn2dem_pos import dn2dem_pos
import dn2dem_pos

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import NullFormatter
import matplotlib
from astropy import units as u

import pathlib
import pickle
import glob

import warnings
warnings.simplefilter('ignore')
matplotlib.rcParams['font.size'] = 16

#Import all of the instrument-specific prep code
import aia_dem_prep
import xrt_dem_prep
import nustar_dem_prep
import visualize_dem_results as vdr


exposure_dict={'Be_thin': [],
                'Be_thick': [],
              'Al_poly': []}

def dodem(time, bl, tr, minT=5.8, maxT=7.5, dT=0.05, 
          xrt=True, aia=True, nustar=True, eis=True,
          plotresp=True, just_prep=False, plotMK=False, name='',
          use_prior_prep=False, default_err=0.2, use_highTprep=False,
          path_to_dodem='./', working_directory='./',
          
          #NuSTAR-related
          nustar_path='./', nuenergies=[[2.5,8]], fpm='A', make_nustar=False, pile_up_corr=False, adjacent_grades=False,
          gtifile='starter_gti.fits', datapath='', regfile='starter_region.reg', edit_regfile=True, use_fit_regfile=False,
          COM_nustar_region=False, compare_fpm=False, combine_fpm=False, nuclobber=False, special_pha='', nuradius=150,
          #XRT-related
          xrt_path='./xrt_for_DEM/', xrt_exclude=[], xrt_factor=2, xrtmethod='First', xrt_exposure_dict=exposure_dict,
          input_xrt_region=[], input_xrt_region_dict=[], real_xrt_err=False, plot_xrt=False,
          #AIA-related
          real_aia_err=False, aia_clobber=False, aia_path='./', aia_exclude=[], aiamethod='Middle', 
          input_aia_region=[], input_aia_region_dict=[], plot_aia=False,
          sunpy_dir='/Users/jessieduncan/sunpy/', 
          errortab='/Users/jessieduncan/ssw/sdo/aia/response/aia_V3_error_table.txt',
          #EIS-related
          eis_path='./', eisbox_bl=[0,0], eisbox_tr=[0,0], contrib_file='chiantipy_gfnt.pickle', 
          edens=1.e+9, eis_exclude=[], fiasco=False, fiascofile='',
          #DEMREG/DEM-related
          mc_in=False, mc_rounds=10, reg_tweak=1.0, max_iter=10, gloci=1, rgt_fact=1.5, 
          dem_norm0=None, nmu=40, emd_int=True, emd_ret=True):
    """
    Wrapper for doing multi-instrument DEMs of a particular region/interval.
    
    Developed from Iain Hannah's DEMREG example notebook:
    https://github.com/ianan/axn_example/blob/main/axn_sep2020_new.ipynb
    
    Also includes an option to output prepped DEM inputs. See later wrapper (run_iterative_wrapper())
    to do a DEM in IDL (via HISSW) using xrt_dem_iterative2.pro.
    
    Initial preparation of data + temperature responses is done using code in the (instrument)_dem_prep modules.
    For documentation on this process for each instrument, see those files. Significant IDL dependence remains,
    some (though not all) implemented through HISSW. 
    

    Keywords - general
    ------------------
    time - interval for DEM (two astropy Time objects)
         FORMAT LIKE, time=[astropy.time.Time('2018-05-29T19:08:00', scale='utc'), 
                                 astropy.time.Time('2018-05-29T19:14:00', scale='utc')]
         
    bl, tr - define regtangular region for DEM (bottom left, top right in arcsec)
        FORMAT LIKE, bl=[-200*u.arcsec,150*u.arcsec]
                      tr=[-50*u.arcsec,300*u.arcsec]
             Used to define XRT, AIA data regions UNLESS input_xrt_region/input_aia_region are set to specify
             region objects. 
                      
    minT, maxT - minimum and maximum temperatures for DEM calculation
                      
    aia, xrt, nustar, eis - set True to include each instrument
                Notes:
                -AIA code is written for use with 6 channels:
                A94, A131, A171, A193, A211, A335
                -NuSTAR code is written to use energy ranges selected via nuenergies keyword
                -XRT can use any number of filters, which is handled automatically. Prepped (level 1) XRT files must
                all be in the directory indicated by xrt_path.
                -eis is under construction
    
    plot - set True to plot XRT, AIA images along the way. 
    
    plotresp - set True to plot temperature responses + loci curves (before DEM calculation).
    
    use_prior_prep - set True to search for an existing output pickle file from a previous DEM run over the same interval,
                        and extract inputs from that file instead of re-running the data prep code.
                        
    default_err - Default error to use (fraction of input). Used for XRT, added in quadrature with NuSTAR statistical 
                    uncertainty, and used for AIA if real_aia_err=False.                    
    
    working_directory - where result images + files will be placed. 
                         Also, location of timestring directory (named like /19-09-00_19-14-00/ after time interval
                         hhmmss_hhmmss). The timestring directory is where prepped data files are placed.
    
    path_to_dodem - location of do-dem installation
    
    Keywords - DEM Calculation
    -------------------------------
    mc_in – set True to do MCMC (repeated DEM calculation with inputs varied within their uncertainty). 
    
    mc_rounds - number of mcmc iterations (intentionally distinct from max_iter, which has meaning within DEMReg)
    
    
        DEMReg Keywords - only enter in DEMReg call
        -------------------------------------------
        (Documentation of these copied from DEMReg documentation.)
    
            gloci - If no dem_norm0 given (or dem_norm0 array of 1s) then set gloci 1 or 0 (default 0) to choose 
                    weighting for the inversion process (L constraint matrix).
                1: uses the min of EM loci curves to weight L.
                0: uses two reg runs - first with L=diag(1/dT) and DEM result from this used to weight L for second run. 

            reg_tweak - The initial normalised chisq target.

            max_iter – The maximum number of iterations to attempt, code iterates if negative DEM is produced. If 
                        max iter is reached before a suitable solution is found then the current solution is returned 
                        instead (which may contain negative values)(Default is only 10 - although non_pos=True will 
                        set as 1).
            rgt_fact – The factor by which rgt_tweak increases each iteration. As the target chisq increases there 
                        is more flexibility allowed on the DEM

            nmu – Number of reg param samples to calculate (default (or <=40) 500 for 0D, 42 for map)

            dem_norm0 – This is an array of length nt which contains an initial guess of the DEM solution providing
                        a weighting for the inversion process (L constraint matrix). The actual values of the 
                        normalisation do not matter, only their relative values. 
                        If no dem_norm0 given then L weighting based on value of gloci (0 is default)

            emd_int – Do the regularization in EMD [cm^-5] instead of DEM [cm^-5 K^-1] space? (default False). In 
                        some circumstances this does seem to help (particularly at higher T), but needs additional 
                        tweaking, so why it is not the default.

            emd_ret – Return EMD solution instead of EMD [cm^-5] instead of DEM [cm^-5 K^-1] (default False)


    Keywords - optional + instrument-specific
    ---------------------------------------
    
    just_prep - for whatever combination of instruments used, return the DEM inputs after they are prepped (don't
                actually do a DEM). 
    
    NuSTAR-related 
    -------------- 
    nuenergies - if using NuSTAR data, provide an energy range (or ranges) to use
            FORMAT LIKE, nuenergies=[2.5,7]
                            OR
                         nuenergies=[[2.5,3.5], [3.5,6.],[6.,10.]]
           
            
    datapath – NuSTAR data directory (location of event_cl, etc. directories). If default downloaded data is used, will
                be named after the OBSID. E.g.: '/Users/jessieduncan/nustar/may-2018/5_29test/80410203001/'            
            
    nustar_path –   Location of timestring directory (named like /19-09-00_19-14-00/ after time interval hhmmss_hhmmss).
                    The timestring directory is where NuSTAR spectral data products, etc. are either found or placed.            
            
    make_nustar - Set True to make nustar spectral data products if they are not found. If you want to use this, you must
                    also provide inputs to gtifile, datapath, regfile (see below).
                          
                   For more documentation, see nustar_dem_prep.py 
                   
    nuclobber - set True to ignore previously-made spectral data products + re-do them all.
                   
    nustar_path –   Location of timestring directory (named like /19-09-00_19-14-00/ after time interval hhmmss_hhmmss).
                    The timestring directory is where NuSTAR spectral data products, etc. are either found or placed.
                    
    gtifile – NuSTAR GTI file – will be edited based on time input. Can find one to use in event_cl directory post 
                initial pipeline (see datapath). 


    special_pha - to enter special PHA files (like to do nonthermal subtraction), set this equal to a path to 
                   their location. See nustar_dem_prep.find_nuproducts() to see expected file naming structure.
                NOTE: each file name must contain the fpm (capitalized) (i.e. 'A' or 'B'), as well as the time interval
                in the format expected ('hh-mm-ss_hh-mm-ss').
                   
                
        Region-Related
        ---------------

        regfile – If edit_regfile=True, this is used as a file template to create an optimized region file. If 
                    edit_regfile=False, this is used as the actual region file for making spectral data products. 

        edit_regfile – Bool: do we want to optimize the region file, or use the one we input?

        use_fit_regfile - If you've previously run the DEM (or something else) to generate a region file via region 
                            fitting, use that file instead of re-running (because it takes a while). Use with
                            edit_regfile=False. (expects a specific name and location; to use a generic user-input 
                            region file set regfile keyword, and edit_regfile=False).
                            
        COM_nustar_region - Set True to use circular region around center of mass of NuSTAR data (ignoring chip gap).
                            Use first with edit_regfile=True (to make new region file). Then afterwards, with
                            edit_regfile=False. 
                            
        nuradius -  radius (in arcseconds) of circular region, for use in the center-of-mass case. Default 150".
        
        Region keywords examples:
            -Use a desired region in an existing file (starter_region.reg):
                    regfile='starter_region.reg', edit_regfile=False
                    
            -Make a new region file, fitting to maximize NuSTAR emission in region (takes a while):
                Note: regfile is just used as a file template.
                    regfile='starter_region.reg', edit_regfile=True, COM_nustar_region=False 
                    
            -Use a region file made from a prior run fit region like the above ^:
                    edit_regfile=False, use_fit_regfile=True
                    
            -Make a new region rule, of 150" radius around NuSTAR center of mass:
                Note: regfile is just used as a file template.
                    regfile='starter_region.reg', edit_regfile=True, COM_nustar_region=True, nuradius=150
                    
            -Use a region file made from a prior run COM region like the above ^:
                Note: regardless of what you set nuradius, the region will have the radius used when last 
                making a COM region file.
                    regfile='starter_region.reg', edit_regfile=False, COM_nustar_region=True
                    

        Pileup/Grade-Related
        ---------------
          
        adjacent_grades - Set True to use grade 0-4 (rather than just grade 0) NuSTAR events 
        
        pile_up_corr - Set True to do a pileup correction. Current method: 
                If adjacent_grades==True: (grades 0-4) - 5/4*(grades 21-24)
                If adjacent_grades==False: (grade 0) - 1/4*(grades 21-24)
                
                
        FPM-Related
        ------------
        
        fpm - if using NuSTAR data, provide a focal plane module ('A' or 'B'). This will allow the correct NuSTAR
            response files to be loaded. (STRING INPUT)
        
        compare_fpm - Make both fpma, fpmb products, use whichever has more NuSTAR emission in its optimal region.
                        (overrides fpm setting)

        combine_fpm - Combine fpma, fpmb data and responses (obviously overrides fpm).
        
    
    xrt-related
    -----------
    xrt_path - if using XRT, path to location of prepped (level 1) XRT files. 
    
    xrt_exclude - list of XRT filters to NOT use.
    
    xrt_factor - multiply all xrt responses by the input factor (default 2).
    
    xrtmethod - Set:
                    - xrtmethod='First' : Use first XRT file found in time interval for each filter
                    - xrtmethod='Average' : Take all files for each filter + average data.
                    
    xrt_exposure_dict - Can enter your own dictionary with a "maximum exposure time" condition for each
                        xrt filter combination in use. Good for cases where some (but not all) files 
                        from a given filter combination are saturated!
                        
    input_xrt_region - Set to use a specific region (via region object) rather than the submap defined by bl, tr.
                            Currently, only set up for "rectangle" or "circle" regions... but other types could 
                            easily be added in xrt_dem_prep.load_xrt_filter().
    
    input_xrt_region_dict - Dictionary of input values for generating region object for each XRT map. 
                                EXAMPLE: (rectangular region)
                                region_data = {'center': (815, -543)*u.arcsec,
                                              'width': 20*u.arcsec,
                                              'height': 60*u.arcsec, 
                                              'angle': 30*u.deg}
                                              
    real_xrt_err - Set True to estimate errors using the expression from Lee et al. 2017, in quadrature with 10% error.
                   Set False to used default_err*xrt input
                                                               
    
    aia-related
    -----------
    
    real_aia_err - set True to use aiapy estimate_error to extract channel-specific estimated 
                    uncertainties (rather than using flat default_err % value). Note you must point to valid aia 
                    error tables (see aia_dem_prep.py).
    
    aia_exclude - wavelength channels to NOT use. Note the default is to use [94, 131, 171, 193, 211, 335].
                    Format should be a list of integers (not strings!). 
                  NOTE: under the assumption that all of these channels will normally be useful in a given NuSTAR/AIA
                  DEM, data for all six will be prepared, and the exclude step comes right before the data is
                  returned for use in the DEM. A typical workflow would be to prep everything first, make an 
                  all-aia DEM, then perhaps remove certain AIA channels if useful in analysis. 
                  
    
    aia_path –   Location of timestring directory (named like /19-09-00_19-14-00/ after time interval hhmmss_hhmmss).
                    The timestring directory is where AIA prepped data files are placed. 
                    
    aiamethod - Set 
                    - aiamethod='Middle' : choose one AIA file near the middle of the time interval to use for 
                                            each wavelength. 
                    – aiamethod='Average' : take AIA file samples every 30s during the interval, and average
                                            to get a DEM-input value. 
                   - aiamethod='Auto' : If the DEM time interval is > 1 minute in length, use "average". If less,
                                           use 'Middle'. Set when doing sets of DEMs over variable time intervals.
                                            
    input_aia_region - Set to use a specific region (via region object) rather than the submap defined by bl, tr.
                            Currently, only set up for "rectangle" or "circle" regions... but other types could 
                            easily be added in aia_dem_prep.map_to_dn_s_px().
    
    input_aia_region_dict - Dictionary of input values for generating region object for each AIA map. 
                                EX: (rectangular region)
                                region_data = {'center': (815, -543)*u.arcsec,
                                              'width': 20*u.arcsec,
                                              'height': 60*u.arcsec, 
                                              'angle': 30*u.deg}
    
    eis-related - UNDER CONSTRUCTION
    -----------
    eis_path - location of eis data files (and fit files, and contribution function file – which are made 
                if not found). 
    eis_exclude - EIS ions to NOT use. Note the default is to use: 
                    ['si_07', 'fe_11', 'fe_09', 'fe_10', 'fe_11', 'fe_12', 'fe_13', 'fe_14', 'fe_16', 'ca_15']
                    Note that there are two fe_11 lines – currently, you can exclude/not exclude by the ion. So
                    if you don't want the first fe_11, you will remove both. 
    
    eisbox_bl, eisbox_tr - bounding box for which section of EIS intensity map to use (spatial selection)
    contrib_file - name of saved file (pickle) which contains EIS-line contribution functions. If file 
                    doesn't exist it's made.
    edens - electron density (for calculating contribution functions)
    
    """
    
    if mc_in==True and just_prep==True:
        print('We will not get to DEMReg MCMC since we are exiting with data products for input elsewhere')
    
    if aia==False and xrt==False:
        print('Not set up to do a DEM without EITHER aia or xrt!')
        return
            
    print('Start Time: ', time[0])
    print('Stop Time: ', time[1])
    
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
    #print(timestring)
    
    dem_path = pathlib.Path(working_directory) / timestring
    if not dem_path.exists():
        dem_path.mkdir()

    #Name the file where we will store the DEM inputs and results. 
    if mc_in:
        picklefile = working_directory+timestring+'/'+timestring+'_'+str(minT)+'_'+str(maxT)+'_'+name+'_MC_DEM_result.pickle'
        if use_highTprep or use_prior_prep:
            gg = glob.glob(working_directory+timestring+'/'+timestring+'_'+str(minT)+'_*_'+name+'_MC_DEM_result.pickle')
            gg.sort(reverse=True)
    else:
        picklefile = working_directory+timestring+'/'+timestring+'_'+str(minT)+'_'+str(maxT)+'_'+name+'_DEM_result.pickle'
        if use_highTprep or use_prior_prep:
            gg=glob.glob(working_directory+timestring+'/'+timestring+'_'+str(minT)+'_'+str(7.2)+'_'+name+'_DEM_result.pickle')
            gg.sort(reverse=True)
        
    if just_prep:
        picklefile=working_directory+\
                    timestring+'/'+timestring+'_'+str(minT)+'_'+str(maxT)+'_'+name+'_iterative_DEM_result.pickle'
        if use_highTprep or use_prior_prep:
            gg = glob.glob(working_directory+\
                           timestring+'/'+timestring+'_'+str(minT)+'_'+str(7.2)+'_'+name+'_iterative_DEM_result.pickle')
            gg.sort(reverse=True)
            
            
    #==========================================================================        
    #========================================================================== 
    #If reusing DEM prep...
    #========================================================================== 
    #========================================================================== 
    
    if use_prior_prep:
        #Using the prior result from this temperature range takes priority over using the default temperature
        #range if both are set.
        use_highTprep=False
    
    prior_run = pathlib.Path(picklefile)
    if use_prior_prep and prior_run.exists() or use_highTprep and bool(gg):
        print('')
        print('Using inputs from prior DEM run, saved in this file: ')
        if use_prior_prep:
            with open(picklefile, 'rb') as f:
                 data = pickle.load(f)
            print(prior_run)
        if use_highTprep:
            with open(gg[0], 'rb') as f:
                 data = pickle.load(f)
            print(gg[0])
            
        #Define DEM inputs, saved in data file.    
        chanax=data['chanax']
        dn_in=data['dn_in']
        edn_in=data['edn_in']
        trmatrix=data['trmatrix']
        temps=data['ts_']
        
        #But editing to use the currently specified instruments (not necessarily just 
        #the ones used in the saved DEM). 
        nuindices = [i for i in range(0, len(chanax)) if 'keV' in chanax[i]]
        xrtindices = [i for i in range(0, len(chanax)) if 'e-' in chanax[i]]
        xrtindices.extend([i for i in range(0, len(chanax)) if 'l-' in chanax[i]]) 
        aiaindices = [i for i in range(0, len(chanax)) if len(chanax[i]) <= 4]
        
        if bool(aia_exclude):
            exclusions = ['A'+str(a) for a in aia_exclude]
            for e in exclusions:
                exindex=[i for i in range(0, len(chanax)) if e == chanax[i]]
                if bool(exindex):
                    aiaindices.remove(exindex[0])
                
        if bool(xrt_exclude):
            for x in xrt_exclude:
                print(chanax)
                exindex=[i for i in range(0, len(chanax)) if x == chanax[i]]
                if bool(exindex):
                    xrtindices.remove(exindex[0])
            
        
        indices=[]
        if aia:
            indices.extend(aiaindices)
        if xrt:
            indices.extend(xrtindices)
        if nustar:
            indices.extend(nuindices)
            
        chanax=[chanax[i] for i in indices]
        dn_in=[dn_in[i] for i in indices]
        edn_in=[edn_in[i] for i in indices]
        trmatrix=trmatrix[:,indices]
      
        nutrs = [c for c in chanax if 'keV' in c]
        nf = len(chanax)
        print('')
        for i in np.arange(len(dn_in)):
            print(chanax[i],':    ',"{0:.2f}".format(dn_in[i]),
                  "  {0:.2f}".format(edn_in[i]),
                  " {0:.0f}".format(100.*edn_in[i]/dn_in[i]),'%')

    
    #========================================================================== 
    #========================================================================== 
    
    #========================================================================== 
    #========================================================================== 
    #If not using prior data prep...
    #========================================================================== 
    #========================================================================== 
    
    
    else:

        #FOR EACH INSTRUMENT, load in: 
        #    -The observed counts (DN/s/px)
        #.   -The list of types of observation (channels, filters, energy ranges, etc) 
        #.   -The temperature response
        #.   -The corresponding list of temperatures


        #data list and instrument name list (to be updated)
        dn_in=[]
        chanax=[]
        #Number of instruments (to be updated)
        nf=0

        if aia:
            #Load in AIA data, uncertainty(if set), labels, temperature response, and corresponding temperatures

            if aiamethod=='Auto':
                dur = (time[1]-time[0]).to(u.s).value
                if dur > 60:
                    aiamethod='Average'
                else:
                    aiamethod='Middle'

            res = aia_dem_prep.load_aia(time, bl, tr, plot=plot_aia, aia_exclude=aia_exclude, aia_path=working_directory, 
                                        method=aiamethod, input_region=input_aia_region, aia_clobber=aia_clobber,
                                        input_aia_region_dict=input_aia_region_dict, real_aia_err=real_aia_err,
                                       sunpy_dir=sunpy_dir,errortab=errortab,path_to_dodem=path_to_dodem)
            if res is None:
                print('Something is wrong; Not using AIA.')
                print('')
                aia=False
            if aia:
                if real_aia_err:
                    aia_dn_s_px, aia_err_dn_s_px, chans, aia_tr, aia_logt = res
                else:
                    aia_dn_s_px, chans, aia_tr, aia_logt = res
                nf+=len(chans)
                for i in range(0, len(chans)):
                    dn_in.append(aia_dn_s_px[i])
                    chanax.append(chans[i])
                    
        if xrt:
            print('')
            #Load in XRT data, uncertainty(if set), labels, temperature response, and corresponding temperatures
            res = xrt_dem_prep.load_xrt(xrt_path, time, bl, tr, xrt_exclude=xrt_exclude, plot=plot_xrt, method=xrtmethod,
                                       exposure_dict=xrt_exposure_dict, input_xrt_region_dict=input_xrt_region_dict,
                                       input_xrt_region=input_xrt_region, real_xrt_err=real_xrt_err,
                                        path_to_dodem=path_to_dodem)
            if res is None:
                print('Something is preventing use of XRT.')
                print('')
                xrt=False
            if xrt:
                if real_xrt_err:
                    xdnspxs, xrt_errs, filters, xrt_tr, xrt_logt  = res
                else:    
                    xdnspxs, filters, xrt_tr, xrt_logt  = res
                if bool(xrt_factor):
                    xrt_tr*=xrt_factor
                    #xrt_tr*=2.5
                xrtstart=nf
                nf+=len(filters)
                xrtstop=nf
                for i in range(0, len(filters)):
                    dn_in.append(xdnspxs[i])
                    if bool(xrt_factor):
                        chanax.append(filters[i]+' x'+str(xrt_factor))
                    else:
                        chanax.append(filters[i])


#         if eis:
#             print('')
#             res = eis_dem_prep.load_chianti_and_eis(eis_path=eis_path, minT=minT, maxT=maxT, dT=dT,
#                                       eisbox_bl=eisbox_bl, eisbox_tr=eisbox_tr, contrib_file=contrib_file, 
#                                       edens=edens, eis_exclude=eis_exclude, fiasco=fiasco,fiascofile=fiascofile)
#             phot_int, linelabels, chi_gfnt, gfnt_temps = res

#             warren_compare=True
#             only_goodlines=False
#             use_warren_int=True

#             if bool(eis_exclude) and warren_compare:
#                 print('')
#                 print('Warren comparison only works well when not using aia_exclude to exclude any ions.')
#                 print('Thus, skipping it.')
#                 warren_compare=False



#             if warren_compare:
#                 #Copied from table 1, Warren (2020)
#                 warren_int = np.array([40.02, 40.76, 27.16, 271.07, 456.01, 993.64, 721.94, 918.21, 866.55, 1566.44,
#                               11018.92, 1321.13, 138.04, 71.48])

#                 warren_lines = np.array([275.368, 188.497, 197.862, 184.536, 188.216, 195.119, 202.044, 
#                        203.826, 270.519, 264.787, 284.160, 262.984, 193.874, 200.972])

#                 #edit_inds = [0,1,2,3,4,5,6,8,11,13]
#                 edit_inds = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]

#                 warren_comp_int = warren_int[edit_inds]
#                 warren_comp_lines = warren_lines[edit_inds]

#                 #phot_warren_comp = eis_dem_prep.convert_intensity(warren_comp_lines, warren_comp_int, ergstart=True)
#                 phot_warren_comp=warren_comp_int

#                 print('Compare intensities:')
#                 print('Ion -- Line -  Warren Value –––––– Data Value -–––––– Ratio')
#                 goodlines=[]
#                 for i in range(0, len(edit_inds)):
#                     ratio = phot_warren_comp[i]/phot_int[i]
#                     print(linelabels[i], warren_comp_lines[i], phot_warren_comp[i], phot_int[i], ratio)
#                     if ratio > 0.5: # and ratio < 5:
#                         #print('Under a half order off: ', warren_comp_lines[i], ratio)
#                         goodlines.append(i)
#                         #print('')
#                 print('')

#                 #If true, use warren-reported intensities for each line rather than what we find
#                 if use_warren_int==True:
#                     #phot_int = phot_warren_comp
#                     phot_int = warren_comp_int

#                 #if true, only use lines fititng the criteria defined above (written to exclude lines
#                 #with intensities much larger or smaller than the corresponding warren value). 
#                 if only_goodlines==True:
#                     phot_int = phot_int[goodlines]
#                     linelabels = [linelabels[i] for i in goodlines]
#                     chi_gfnt = chi_gfnt[goodlines, :]


#             if res is None:
#                 print('Something is wrong; Not using EIS.')
#                 print('')
#                 eis=False
#             if eis:
#                 print('Doin EIS')
#                 nf+=len(linelabels)
#                 for i in range(0, len(linelabels)):
#                     dn_in.append(phot_int[i])
#                     chanax.append(linelabels[i])


        #Define error in each observed value - default flat %
        edn_in=list(default_err*np.copy(dn_in))

        #If set, switch to real aia errors
        if real_aia_err and aia:
            edn_in[0:len(aia_err_dn_s_px)] = aia_err_dn_s_px
            
        #if set, switch to real xrt errors
        if real_xrt_err and xrt:
            edn_in[xrtstart:xrtstop] = xrt_errs


        if nustar:
            #Load in NuSTAR data, uncertainty, labels, temperature response, corresponding temperatures, and fpm
            print('')

            if combine_fpm:
                res = nustar_dem_prep.combine_fpm(time, nuenergies, working_directory, make_nustar=make_nustar,
                                              gtifile=gtifile, datapath=datapath, regfile=regfile, 
                                              edit_regfile=edit_regfile, use_fit_regfile=use_fit_regfile,
                                                 nofit=COM_nustar_region, clobber=nuclobber, pile_up_corr=pile_up_corr,
                                                  default_err=default_err, special_pha=special_pha,
                                                 adjacent_grades=adjacent_grades, nuradius=nuradius,
                                                  path_to_dodem=path_to_dodem)
            else:

                res = nustar_dem_prep.load_nustar(time, nuenergies, working_directory, fpm, make_nustar=make_nustar,
                                                  gtifile=gtifile, datapath=datapath, regfile=regfile, 
                                                  edit_regfile=edit_regfile, compare_fpm=compare_fpm,
                                                  use_fit_regfile=use_fit_regfile, pile_up_corr=pile_up_corr,
                                                 nofit=COM_nustar_region, clobber=nuclobber, 
                                                 default_err=default_err, special_pha=special_pha,
                                                 adjacent_grades=adjacent_grades, nuradius=nuradius, 
                                                 path_to_dodem=path_to_dodem)
            if res is None:
                print('Something is wrong; Not using NuSTAR.')
                print('')
                nustar=False

            if nustar:
                rate, erate, nutrs, nu_tresp, nu_logt, fpm = res
                nf+=len(nutrs)
                for i in range(0, len(nutrs)):
                    dn_in.append(rate[i])
                    edn_in.append(erate[i])
                    if combine_fpm==False:
                        chanax.append(nutrs[i]+'-'+fpm)
                    else:
                        chanax.append(nutrs[i])

        #======================================================

        #======================================================
        #Print rates and errors

        print('')
        for i in np.arange(len(dn_in)):
            print(chanax[i],':    ',"{0:.2f}".format(dn_in[i]),
                  "  {0:.2f}".format(edn_in[i]),
                  " {0:.0f}".format(100.*edn_in[i]/dn_in[i]),'%')


        #======================================================

        #======================================================
        #Define the temperatures which will be used as DEM input 
        #+ make unified response matrix:


        if aia: 
            #Temperatures to use: AIA response temperatures
            temps = aia_logt
            nt=len(temps)
            trmatrix=np.zeros((nt,nf))
            #Add aia responses to response matrix
            for i in range(0,len(chans)):
                trmatrix[:,i]=aia_tr[i] #units counts*cm^5/s
            added=len(chans)
            #If including, add XRT responses to response matrix
            if xrt==True:
                for f in range(0, len(filters)):
                    #Interpolate XRT response at input temperatures
                    xrtint=10**np.interp(temps,xrt_logt,np.log10(xrt_tr[f]))
                    trmatrix[:,(added+f)]=xrtint
                added+=len(filters)
            #If including, add NuSTAR responses to response matrix
            if nustar==True:
                #Interpolate NuSTAR response at input temperatures
                for i in range(0, len(nutrs)):
                    intp_ntresp=10**np.interp(temps,nu_logt,np.log10(nu_tresp[:,i]))
                    trmatrix[:,(added+i)]=intp_ntresp
                added+=len(nutrs)
#             if eis==True:
#                 #Interpolate EIS G(n,T) at input temperatures
#                 for i in range(0, len(linelabels)):
#                     intp_egfnt=np.interp(temps,np.log10(gfnt_temps),chi_gfnt[i,:])
#                     trmatrix[:,(added+i)]=intp_egfnt



        if aia==False and xrt==True:
            added=0
            #Temperatures to use: XRT response temperatures
            temps = xrt_logt
            nt=len(temps)
            trmatrix=np.zeros((nt,nf))
            for f in range(0, len(filters)):
                trmatrix[:,f]=xrt_tr[f]
            added+=len(filters)
            #If including, add NuSTAR responses to response matrix
            if nustar==True:
                #Interpolate NuSTAR response at input temperatures
                for i in range(0, len(nutrs)):
                    intp_ntresp=10**np.interp(temps,nu_logt,np.log10(nu_tresp[:,i]))
                    trmatrix[:,(added+i)]=intp_ntresp
                added+=len(nutrs)
#             if eis==True:
#                 #Interpolate EIS G(n,T) at input temperatures
#                 for i in range(0, len(linelabels)):
#                     intp_egfnt=np.interp(temps,np.log10(gfnt_temps),chi_gfnt[i,:])
#                     trmatrix[:,(added+i)]=intp_egfnt

        #=================================================================================================
        #=================================================================================================
        #=================================================================================================
    
    
    #======================================================        
    # Set up some colours for plotting
    if nustar:
        clrs = make_DEM_colors(aia=aia, aia_exclude=aia_exclude, xrt=xrt, xrt_exclude=xrt_exclude,
                                        nustar=nustar, nustar_number=len(nutrs),
                                        eis=eis, eis_exclude=eis_exclude)
    else:
        clrs = make_DEM_colors(aia=aia, aia_exclude=aia_exclude, xrt=xrt, xrt_exclude=xrt_exclude,
                                        nustar=nustar, eis=eis, eis_exclude=eis_exclude)
            

    #======================================================
    
    
    #======================================================
    #EXTRACT PREPPED DATA IF SET
    
    if just_prep==True:
        return dn_in, edn_in, trmatrix, temps, chanax
    
    #======================================================
    
    #======================================================
    #Make some response plots if set
    
    if plotresp==True:
        
        if eis==False:
            #Doesn't make sense to do with EIS, since "response" units different
    
            # Plot all the temperature responses
            fig = plt.figure(figsize=(9, 7))
            plt.rcParams.update({'font.size': 16,'font.family':"sans-serif",\
                                 'font.sans-serif':"Arial",'mathtext.default':"regular"})
            for i in np.arange(len(chanax)):
                plt.semilogy(temps,trmatrix[:,i],label=chanax[i],color=clrs[i],lw=4)
            plt.xlabel('$\mathrm{\log_{10}T\;[K]}$')
            plt.ylabel('$\mathrm{Response\;[(DN\;px^{-1}\;\;or\;\;Count)\;s^{-1}\;cm^5]}$')
            #plt.ylim([2e-29,5e-24])
            plt.xlim([5.6,7.2])
            plt.title('Temperature Response of Instruments Used for DEM')
            plt.legend(ncol=2,prop={'size': 14})
            plt.rcParams.update({'font.size': 16})
            plt.grid(True,which='both',lw=0.5,color='gainsboro')
            plt.show()
    
    
        # Plot the EM Loci
        print('EM Loci are just the input rate divided by the temperature response:')

        plt.rcParams.update({'font.size': 16,'font.family':"sans-serif",\
                             'font.sans-serif':"Arial",'mathtext.default':"regular"})
        fig = plt.figure(figsize=(9, 7))
        for i in np.arange(len(chanax)):
            plt.semilogy(temps,dn_in[i]/trmatrix[:,i],label=chanax[i],color=clrs[i],lw=4)
        plt.ylim([1e23,1e30])
        plt.xlim([5.8,7.2])
        plt.xlabel('$\mathrm{\log_{10}T\;[K]}$')
        plt.ylabel('$Emission\;Measure\;[cm^{-5}]$')
        plt.legend(ncol=2,prop={'size': 10})
        plt.grid(True,which='both',lw=0.5,color='gainsboro')
        fig.show()
    
    #======================================================
    
    #======================================================
    #Make DEM
    
    #Temperature bins
    dlogt=dT
    #Array of temperatures (not in log space)
    temps70=10**np.arange(minT,maxT+dlogt,dlogt)
    #For each temperature, take the mean of that temperature and the next one up (and keep all in an array, log space)
    mlogt70=([np.mean([(np.log10(temps70[i])),np.log10((temps70[i+1]))]) \
            for i in np.arange(0,len(temps70)-1)])
    


#     print('')
#     print('DEM Inputs:')
#     print('-Rates, Rate errors')
#     print('-Temperature response matrix, temperature response temperatures')
#     print('-Input temperatures')

#     print('')
#     print('Outputs:')
#     print('-DEM, vertical errors in DEM')
#     print('-Horizontal errors on temperature')
#     print('-Chisq')
#     print('-Simulated counts in each channel, based on DEM results + response')
    
    #DO DEM
    dem70o,edem70o,elogt70o,chisq70o,dn_reg70o\
        =dn2dem_pos.dn2dem_pos(np.array(dn_in), np.array(edn_in), trmatrix, temps, temps70, gloci=gloci, emd_int=emd_int,
                               emd_ret=emd_ret, reg_tweak=reg_tweak, max_iter=max_iter, rgt_fact=rgt_fact, 
                               dem_norm0=dem_norm0, nmu=nmu)
    
    if np.any((dem70o < 0)):
        print('')
        print('WARNING: DEM SOLUTION HAS NEGATIVE-DEM TEMP. BINS!!!')
        print('Not saving!! Re-TRYING')
        tries=5
        while tries > 0:
            dem70o,edem70o,elogt70o,chisq70o,dn_reg70o\
                    =dn2dem_pos.dn2dem_pos(np.array(dn_in), np.array(edn_in), trmatrix, temps, temps70, gloci=gloci,
                                           emd_int=emd_int, emd_ret=emd_ret, reg_tweak=reg_tweak, max_iter=max_iter,
                                           rgt_fact=rgt_fact, dem_norm0=dem_norm0, nmu=nmu)
            if np.any((dem70o < 0)) == False:
                print('Now the solution has no <0 bins, continuing.')
                continue
            tries-=1

        if np.any((dem70o < 0)):
            print('After five tries we still have no non-zero solution, quitting.')
            return
                
            print('')
            

            
    #======================================================
    
    #======================================================
    #Get temperatures sorted for when we plot
    
    if plotMK == True:
        #T in MK
        mkt = 10**np.array(mlogt70)/1e6
        mkt_ = 10**np.array(temps)/1e6
    else:
        mkt = np.array(mlogt70)
        mkt_ = np.array(temps)
        
    #======================================================
      
    #Save inputs (outputs added later)
    data = {'time_interval': time,
            'nuenergies': nuenergies,
            'nufpm': fpm,
            'temp_interval': [minT, maxT],
            'plotMK': plotMK,
            'ts': mkt,
            'ts_': mkt_,
            'trmatrix': trmatrix,
            'dn_in': dn_in,
            'edn_in': edn_in,
            'chanax': chanax
               }

    #======================================================
    #Do MCMC if set (& plot, & finish)
    
#     fig = plt.figure(figsize=(9, 12), tight_layout = {'pad': 1})
    
    if mc_in==True:
        print('')
        print('Doing '+str(mc_rounds)+' iterations of DEMReg with input varied within uncertainty!')
        
        dnins = np.zeros((mc_rounds, len(dn_in)))
        dnins[:] = np.nan
        demouts = np.zeros((mc_rounds, len(dem70o)))
        demouts[:] = np.nan
        chisqs = np.zeros(mc_rounds)
        chisqs[:] = np.nan
        dnouts = np.zeros((mc_rounds, len(dn_reg70o)))
        dnouts[:] = np.nan
        
        while mc_rounds > 0:
            
            error=(np.random.rand(len(dn_in))*2-1)*edn_in
            dn_in_now = dn_in+error
            #DO DEM
            dem70oMC,edem70oMC,elogt70oMC,chisq70oMC,dn_reg70oMC\
                =dn2dem_pos.dn2dem_pos(np.array(dn_in_now), np.array(edn_in), trmatrix, temps, temps70, gloci=gloci,
                                       emd_int=emd_int,
                                       emd_ret=emd_ret, reg_tweak=reg_tweak, max_iter=max_iter, rgt_fact=rgt_fact, 
                                       dem_norm0=dem_norm0, nmu=nmu)
            
            if np.any((dem70oMC <= 0)):
                print('')
                print('WARNING: ITERATION ', mc_rounds, ' HAS 0 or NEGATIVE-DEM TEMP. BINS!!!')
                print('Not saving this iteration!!')
                print('')
                mc_rounds-=1
                continue
                
#             plt.semilogy(mkt,dem70oMC)
            
            dnins[mc_rounds-1, :] = dn_in_now
            demouts[mc_rounds-1, :] = dem70oMC
            dnouts[mc_rounds-1, :] = dn_reg70oMC
            chisqs[mc_rounds-1] = chisq70oMC
            mc_rounds-=1
                           
        #Mean, min, max value of all output DEM solutions
        #Mask out nans (never updated rows, from DEMs that failed to find a positive solution)
        demouts_ma = np.ma.masked_invalid(demouts)
        min_outs = demouts_ma.min(axis=0) #np.amin(demouts, axis=0)
        max_outs = demouts_ma.max(axis=0) #np.amax(demouts, axis=0)
        
        #Mean, min, max of all MC DN outputs
        #Mask out nans (never updated rows, from DEMs that failed to find a positive solution)
        dnouts_ma = np.ma.masked_invalid(dnouts)
        min_dnouts = dnouts_ma.min(axis=0) #np.min(dnouts, axis=0)
        max_dnouts = dnouts_ma.max(axis=0) #np.max(dnouts, axis=0)
        
        
        asymmetric_error = [(dem70o-min_outs), (max_outs-dem70o)]
        dem = dem70o
        asymmetric_error_dn = [(dn_reg70o-min_dnouts)/dn_in, (max_dnouts-dn_reg70o)/dn_in]
        dn_reg=dn_reg70o
        edn_string = 'Uncertainties from MCMC output'
        
        
        outputs = {'DEM': dem,
            'dn_reg': dn_reg,
            'edem': asymmetric_error,
            'edn': asymmetric_error_dn, 
            'chisq': chisq70o,
            'edn_string': edn_string,
            'fill_color': 'lightcoral' 
            }
        
        data = data | outputs
        
        #Save to file defined at begining of function
        with open(picklefile, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
            
        vdr.plot_DEM(data, fill_color='lightcoral',
                              title=working_directory+timestring+'/'+timestring+'_'+str(minT)+str(maxT)+name+'_MC')
        
        
        return picklefile

    
    #======================================================
    
    #======================================================
    #SAVE RESULTS + MAKE PLOT - if not doing MCMC
    
    
    asymmetric_error_dn = [(dn_reg70o/dn_in)-dn_reg70o/(np.array(dn_in)+np.array(edn_in)),
                        dn_reg70o/(np.array(dn_in)-np.array(edn_in))-(dn_reg70o/dn_in)]
                        
    edn_string = 'Uncertainties from data input'
    
    outputs = {'DEM': dem70o,
            'dn_reg': dn_reg70o,
            'edem': edem70o,
            'edn': asymmetric_error_dn, 
            'chisq': chisq70o,
            'edn_string': edn_string,
            'xdem_error': elogt70o,
            'fill_color': 'moccasin'   
            }
    
    data = data | outputs

    #Save to file defined at begining of function
    with open(picklefile, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

    vdr.plot_DEM(data, fill_color='moccasin',
                          title=working_directory+timestring+'/'+timestring+'_'+str(minT)+'_'+str(maxT)+name)
    return picklefile



        
def make_DEM_colors(aia=True, aia_exclude=[], xrt=True, xrt_exclude=[], nustar=True, nustar_number=3, 
                    eis=True, eis_exclude=[]):
    """
    Makes a list of colors depending on which instruments are being used (so color-channel correspondence remains 
    consistent across DEMs using different instruemnt combinations.
    
    Color table assumes set numbers of inputs per included instrument:
    
    six AIA channels
    two XRT filters
    three NuSTAR energy ranges
    ten EIS lines.
    
    """
    #Highlight EIS:
    aiaclrs=['dodgerblue','dodgerblue','dodgerblue','dodgerblue','dodgerblue','dodgerblue']
    xrtclrs=['gold', 'gold', 'gold']
    eisclrs=['darkgreen','darkcyan','sienna','indianred','darkorange','seagreen','springgreen','mediumaquamarine',
            'turquoise','lawngreen', 'cadetblue','slategray', 'darkslateblue', 'black']
    
    
    #normal
    aiaclrs=['darkgreen','darkcyan','gold','sienna','indianred','darkorange']
    xrtclrs=['darkslateblue', 'dodgerblue', 'cornflowerblue']
    nuclrs=['purple', 'mediumpurple', 'plum', 'violet']
    eisclrs=['seagreen', 'mediumseagreen', 'springgreen', 'green', 'mediumspringgreen', 'mediumaquamarine', 
          'aquamarine', 'turquoise', 'lightseagreen', 'mediumturquoise', 'lawngreen', 'cadetblue', 
             'slategray', 'darkslateblue']
    clrs=[]
    
    if aia==True:
        index=len(aiaclrs)-len(aia_exclude)
        clrs.extend(aiaclrs[0:index])
        
    if xrt==True:
        index=len(xrtclrs)-len(xrt_exclude)
        clrs.extend(xrtclrs[0:index])
        
    if nustar==True:
        clrs.extend(nuclrs)
        
    if eis==True:
        index=len(eisclrs)-len(eis_exclude)
        clrs.extend(eisclrs[0:index])
        
    return clrs





def run_iterative_wrapper(time, bl, tr, minT, maxT, xrt=True, aia=True, nustar=True, eis=False,
                          name='', return_inputs=False,save_inputs=False, plotMK=True,
                          use_prior_prep=False, special_pha='', mc_iter=100, dT=0.051, chi_thresh=0.95, 
                          default_err=0.2, use_highTprep=False,
                          
                          nuenergies=[2.5,7], pile_up_corr=False, adjacent_grades=False,
                          combine_fpm=False, nustar_path='./', COM_nustar_region=False, nuclobber=False,
                          edit_regfile=False, use_fit_regfile=False,
                          
                          aiamethod='Auto', aia_exclude=[], input_aia_region=[], input_aia_region_dict=[], 
                          real_aia_err=False, aia_clobber=False, 
                          
                          xrtmethod='Average', xrt_exposure_dict=[], xrt_factor=2,
                          real_xrt_err=False, input_xrt_region=[],input_xrt_region_dict=[]):

    """
    IDL wrapper - loads in the needed data and responses, then runs IDL procedure (using HISSW) 
    which converts the responses into the desired structure format and calculates the DEM with 
    xrt_dem_iterative2.pro. Finally, makes a nice plot. 
    
    If return_inputs is set True, instead just makes the inputs for a DEM calculation, and returns them
    for use. 
    
    If save_inputs is set True, makes a lot of .csv files with inputs for a DEM calculation.

    """
    import hissw
    
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
    
    #Run dodem with "just_prep=True", in order to prepare data+responses for DEM
    dn_in, edn_in, trmatrix, temps, chanax = dodem(time, bl, tr, minT=minT, maxT=maxT, xrt=xrt, aia=aia, nustar=nustar, 
                                                   eis=eis, plot=False, just_prep=True, pile_up_corr=pile_up_corr,
                                                   adjacent_grades=adjacent_grades,
                                                   nuenergies=nuenergies, #chu=chu, obsid=obsid, fpm=fpm, 
                                                   nustar_path=nustar_path, COM_nustar_region=COM_nustar_region,
                                                   nuclobber=nuclobber, name=name, special_pha=special_pha,
                                                   #xrt_path=xrt_path, xrt_exclude=xrt_exclude,
                                                   combine_fpm=combine_fpm, aiamethod=aiamethod, xrtmethod=xrtmethod,
                                                   xrt_exposure_dict=xrt_exposure_dict,
                                                   xrt_factor=xrt_factor, real_xrt_err=real_xrt_err,
                                                   aia_exclude=aia_exclude,
                                                   edit_regfile=edit_regfile, use_fit_regfile=use_fit_regfile,
                                                  input_xrt_region=input_xrt_region,
                                                  input_xrt_region_dict=input_xrt_region_dict,
                                                  input_aia_region=input_aia_region, aia_clobber=aia_clobber,
                                                  input_aia_region_dict=input_aia_region_dict,
                                                  use_prior_prep=use_prior_prep,
                                                  real_aia_err=real_aia_err, default_err=default_err,
                                                  use_highTprep=use_highTprep)
    
    
    if save_inputs==True:
        #Saving 1D components
        #format='%.6f', 
        np.savetxt(timestring+'_dn_in.csv', 
                   dn_in, 
                   delimiter =", ", 
                   fmt ='% s')
        np.savetxt(timestring+'_edn_in.csv', 
                   edn_in, 
                   delimiter =", ", 
                   fmt ='% s')

        np.savetxt(timestring+'_temps.csv', 
                   temps, 
                   delimiter =", ", 
                   fmt ='% s')

        #Saving 2D response matrix
        with open(timestring+'_trmatrix.txt', 'w') as outfile:  
            json.dump(trmatrix.tolist(), outfile)

        # using the savetxt 
        # from the numpy module
        # to save the channel names list
        np.savetxt(timestring+'_chanax.csv', 
                   chanax,
                   delimiter =", ", 
                   fmt ='% s')

    if return_inputs==True:
        return dn_in, edn_in, trmatrix, temps, chanax
    
    #Do DEM Calculation in IDL using HISSW:
    inputs = {'obs_val': list(dn_in), 'obs_err': list(edn_in), 'temps': list(temps), 
         'obs_index': list(chanax), 'trmatrix': trmatrix.tolist(), 'mc_iter': mc_iter, 'dT': dT, 
          'Name': [timestring], 'min_T': minT, 'max_T': maxT}
    ssw = hissw.Environment(ssw_packages=['hinode/xrt','sdo/aia'], ssw_paths=['xrt', 'aia'])
    ssw_resp = ssw.run('/Users/jessieduncan/dems/calc_iterative.pro', args=inputs)

    #Fetch DEM results + other helpful values from hissw run
    dem_out=ssw_resp['dem_out']
    chisq=ssw_resp['chisq']
    dlogTfac=ssw_resp['dlogtfac']
    mod_obs=ssw_resp['mod_obs']
    og_dem = dem_out[0,:]*dlogTfac
    logt_out = ssw_resp['logt_out']
    
    data = {'time_interval': time,
            'nuenergies': nuenergies,
            'dem_out': ssw_resp['dem_out'],
            'temp_interval': [minT, maxT],
            'plotMK': plotMK,
            'chisq': ssw_resp['chisq'],
            'dlogTfac': ssw_resp['dlogtfac'],
            'mod_obs': ssw_resp['mod_obs'],
            'ts': ssw_resp['logt_out'], 
            'ts_': temps,
            'trmatrix': trmatrix,
            'dn_in': dn_in,
            'edn_in': edn_in,
            'chanax': chanax,
            'mc_iter': mc_iter
            }
    
    data = read_iterative_outputs(data, chi_thresh=chi_thresh, name=name)
    
    temp_interval=data['temp_interval']
    minT=temp_interval[0]
    maxT=temp_interval[1]
    
    picklefile='./'+timestring+'/'+timestring+'_'+str(minT)+'_'+str(maxT)+'_'+name+'_iterative_DEM_result.pickle'
    with open(picklefile, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)   
        
    return picklefile
    

        

def read_iterative_outputs(data, chi_thresh=0.95, name=''):
    
    import scipy
        
    time=data['time_interval']
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
        
        
    og_dem = data['dem_out'][0,:]*data['dlogTfac']
    og_chisq = data['chisq'][0]
    dn_reg = data['mod_obs'][0]
    
    #MCMC Error bars: set chi squared criteria based on DOF (number of measurements - 1, for number of 
    #splines used in DEM calculation - see xrt_dem_iterative2.pro documentation).
    
    #Error bars are minima/maxima of DEM solution values at each temperature across all MCMC runs that
    #satisfy the chisq condition.
    
    csq = scipy.stats.chi2.ppf(chi_thresh, len(data['chanax'])-1)
    print('Chisq Cutoff: ', csq)
    good_chis = np.where(data['chisq'] < csq)[0]
    print('Good Chisq Fraction: ', len(good_chis)/len(data['chisq']))

    lower_error=np.zeros(len(og_dem))
    upper_error=np.zeros(len(og_dem))
    for i in range(0, len(og_dem)):
        all_res = data['dem_out'][:,i]*data['dlogTfac'][i]
        all_res = all_res[good_chis]
        lower_error[i]=og_dem[i]-min(all_res)
        upper_error[i]=max(all_res)-og_dem[i]
    asymmetric_error = [lower_error, upper_error]

    
    
    data['edem'] = asymmetric_error
    data['DEM'] = og_dem
    data['dn_reg'] = dn_reg
    data['fill_color'] = 'lightgreen'
    
    #Minimum and maximum possible residual values for each instrument, based on input data uncertainties.
    asymmetric_error_dn = [(dn_reg/data['dn_in'])-dn_reg/(np.array(data['dn_in'])+np.array(data['edn_in'])),
                        dn_reg/(np.array(data['dn_in'])-np.array(data['edn_in']))-(dn_reg/data['dn_in'])]
    
    
    data['edn'] = asymmetric_error_dn
    data['edn_string'] = 'Uncertainties from data input'
    
    temp_interval=data['temp_interval']
    minT=temp_interval[0]
    maxT=temp_interval[1]
    
    picklefile='./'+timestring+'/'+timestring+'_'+str(minT)+'_'+str(maxT)+'_'+name+'_iterative_DEM_result.pickle'
    
    with open(picklefile, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
        
    vdr.plot_DEM(data, fill_color='lightgreen',
                          title='./'+timestring+'/'+timestring+'_'+str(temp_interval[0])+str(temp_interval[1])+'_iterative')
        
    return data    
    
def high_temp_analysis(time, bl, tr, exposure_dict, datapath, nuenergies, highT=7.2, dT=0.05,
                       gtifile='starter_gti.fits', regfile='starter_region.reg', name2='',
                       xrt=True, aia=True, nustar=True, edit_regfile=False, use_fit_regfile=False,
                       COM_nustar_region=False, nuclobber=False, special_pha='', pile_up_corr=False,
                       adjacent_grades=False,
                       plotMK=False, plot=False, aia_exclude=[], aia_clobber=False,
                       xrt_exclude=[], xrtmethod='Average', xrt_path='./xrt_for_DEM/', xrt_factor=2,
                       input_xrt_region=[], input_xrt_region_dict=[], real_xrt_err=False,
                       aiamethod='Average', input_aia_region=[], input_aia_region_dict=[],
                       demmethod='DEMREG', use_prior_prep=False,
                       real_aia_err=False, default_err=0.2,
                       reg_tweak=1, max_iter=30, rgt_fact=1.5):
   
    
    """
    Do the DEM for a given time interval, etc. four times with different temperature ranges; do 
    some preliminary comparison.
    
    For more detailed analysis of results in cases where we do this over a range of time intervals,
    see compare_dems.DEM_param_evolution().
    
    """                   
    fpm='A+B'  
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    print(timestring, stopstring)
    timestring=timestring+'_'+stopstring   
    
    if COM_nustar_region:
        name='nofitreg'+name2
    else:
        name=name2
        
    #=========================================================================================================
    minT1=5.6
    #maxT1=7.1
    maxT1=highT
        
        
    if demmethod=='DEMREG':                   
        res1 = dodem(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, eis=False, dT=dT,
                    name=name, compare_fpm=False, combine_fpm=True,
                    nuenergies=nuenergies, nustar_path='./', pile_up_corr=pile_up_corr,
                    make_nustar=True, gtifile=gtifile, datapath=datapath, adjacent_grades=adjacent_grades,
                    regfile=regfile, edit_regfile=edit_regfile, use_fit_regfile=use_fit_regfile,
                    COM_nustar_region=COM_nustar_region, nuclobber=nuclobber, special_pha=special_pha,
                    xrt_path=xrt_path, plot=plot, xrt_exclude=xrt_exclude, 
                    xrt_factor=xrt_factor, xrtmethod=xrtmethod, xrt_exposure_dict=exposure_dict,
                    input_xrt_region=input_xrt_region, input_xrt_region_dict=input_xrt_region_dict,
                    real_xrt_err=real_xrt_err,
                    aiamethod=aiamethod, input_aia_region=input_aia_region, aia_clobber=aia_clobber, 
                    input_aia_region_dict=input_aia_region_dict,
                    just_prep=False,
                    gloci=1, mc_in=True, mc_rounds=100, 
                    plotresp=False, reg_tweak=reg_tweak, rgt_fact=rgt_fact, max_iter=max_iter,
                    minT=minT1, maxT=maxT1, plotMK=plotMK, use_prior_prep=use_prior_prep,
                    real_aia_err=real_aia_err, default_err=default_err)
        
    if demmethod=='XRT_ITER':
        res1 = run_iterative_wrapper(time, bl, tr, minT1, maxT1, xrt=xrt, aia=aia, nustar=nustar, eis=False,
                          dT=dT, nuenergies=nuenergies, nustar_path='./', pile_up_corr=pile_up_corr,
                          edit_regfile=edit_regfile, use_fit_regfile=use_fit_regfile, adjacent_grades=adjacent_grades,
                          COM_nustar_region=COM_nustar_region, nuclobber=nuclobber, special_pha=special_pha,
                          xrt_factor=xrt_factor, real_xrt_err=real_xrt_err,
                          aia_exclude=aia_exclude, name=name,
                         plotMK=plotMK, mc_iter=100,
                         combine_fpm=True, aiamethod=aiamethod, input_aia_region=input_aia_region,
                         input_aia_region_dict=input_aia_region_dict, xrtmethod=xrtmethod, aia_clobber=aia_clobber, 
                         input_xrt_region=input_xrt_region, input_xrt_region_dict=input_xrt_region_dict,
                         xrt_exposure_dict=exposure_dict, use_prior_prep=use_prior_prep,
                        real_aia_err=real_aia_err, default_err=default_err)
        
    #========================================================================================================= 
    minT2=5.6
    maxT2=7.0
    
    if demmethod=='DEMREG':
        res2 = dodem(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, eis=False, dT=dT,
                    name=name, compare_fpm=False, combine_fpm=True,
                    nuenergies=nuenergies, nustar_path='./', pile_up_corr=pile_up_corr,
                    make_nustar=True, gtifile=gtifile, datapath=datapath, adjacent_grades=adjacent_grades,
                    regfile=regfile, edit_regfile=edit_regfile,  use_fit_regfile=use_fit_regfile,
                    COM_nustar_region=COM_nustar_region, nuclobber=False, special_pha=special_pha,
                    xrt_path=xrt_path, plot=plot, xrt_exclude=xrt_exclude, 
                    xrt_factor=xrt_factor, xrtmethod=xrtmethod, xrt_exposure_dict=exposure_dict,
                    input_xrt_region=input_xrt_region, input_xrt_region_dict=input_xrt_region_dict,
                    real_xrt_err=real_xrt_err,
                    aiamethod=aiamethod, input_aia_region=input_aia_region, aia_clobber=aia_clobber, 
                    input_aia_region_dict=input_aia_region_dict,
                    just_prep=False,
                    gloci=1, mc_in=True, mc_rounds=100, 
                    plotresp=False, reg_tweak=reg_tweak, rgt_fact=rgt_fact, max_iter=max_iter, 
                    minT=minT2, maxT=maxT2, plotMK=plotMK,use_prior_prep=use_prior_prep,
                    real_aia_err=real_aia_err, default_err=default_err, use_highTprep=True)  
    
    if demmethod=='XRT_ITER':
        res2 = run_iterative_wrapper(time, bl, tr, minT2, maxT2, xrt=xrt, aia=aia, nustar=nustar, eis=False,
                          dT=dT, nuenergies=nuenergies, nustar_path='./', pile_up_corr=pile_up_corr,
                          edit_regfile=edit_regfile, use_fit_regfile=use_fit_regfile, adjacent_grades=adjacent_grades,
                          COM_nustar_region=COM_nustar_region, nuclobber=False, special_pha=special_pha,
                          xrt_factor=xrt_factor, real_xrt_err=real_xrt_err,
                          aia_exclude=aia_exclude, name=name,
                         plotMK=plotMK, mc_iter=100,
                         combine_fpm=True, aiamethod=aiamethod, input_aia_region=input_aia_region,
                         input_aia_region_dict=input_aia_region_dict, xrtmethod=xrtmethod, 
                         input_xrt_region=input_xrt_region, input_xrt_region_dict=input_xrt_region_dict,
                         xrt_exposure_dict=exposure_dict, use_prior_prep=use_prior_prep, aia_clobber=aia_clobber,
                            real_aia_err=real_aia_err, default_err=default_err, use_highTprep=True)

    #=========================================================================================================                 
    minT3=5.6
    maxT3=6.84
    
    if demmethod=='DEMREG':
        res3 = dodem(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, eis=False, dT=dT,
                    name=name, compare_fpm=False, combine_fpm=True,
                    nuenergies=nuenergies, nustar_path='./', pile_up_corr=pile_up_corr,
                    make_nustar=True, gtifile=gtifile, datapath=datapath, adjacent_grades=adjacent_grades,
                    regfile=regfile, edit_regfile=edit_regfile,  use_fit_regfile=use_fit_regfile,
                    COM_nustar_region=COM_nustar_region, nuclobber=False, special_pha=special_pha,
                    xrt_path=xrt_path, plot=plot, xrt_exclude=xrt_exclude, 
                    xrt_factor=xrt_factor, xrtmethod=xrtmethod, xrt_exposure_dict=exposure_dict,
                    input_xrt_region=input_xrt_region, input_xrt_region_dict=input_xrt_region_dict,
                    real_xrt_err=real_xrt_err,
                    aiamethod=aiamethod, input_aia_region=input_aia_region, aia_clobber=aia_clobber,
                    input_aia_region_dict=input_aia_region_dict, 
                    just_prep=False,
                    gloci=1, mc_in=True, mc_rounds=100, 
                    plotresp=False, reg_tweak=reg_tweak, rgt_fact=rgt_fact, max_iter=max_iter, 
                    minT=minT3, maxT=maxT3, plotMK=plotMK, use_prior_prep=use_prior_prep,
                    real_aia_err=real_aia_err, default_err=default_err, use_highTprep=True) 
    
    if demmethod=='XRT_ITER':
        res3 = run_iterative_wrapper(time, bl, tr, minT3, maxT3, xrt=xrt, aia=aia, nustar=nustar, eis=False,
                          dT=dT, nuenergies=nuenergies, nustar_path='./', pile_up_corr=pile_up_corr,
                          edit_regfile=edit_regfile, use_fit_regfile=use_fit_regfile, adjacent_grades=adjacent_grades,
                          COM_nustar_region=COM_nustar_region, nuclobber=False, special_pha=special_pha,
                          xrt_factor=xrt_factor, real_xrt_err=real_xrt_err, 
                          aia_exclude=aia_exclude, name=name,
                         plotMK=plotMK, mc_iter=100,
                         combine_fpm=True, aiamethod=aiamethod, input_aia_region=input_aia_region,
                         input_aia_region_dict=input_aia_region_dict, xrtmethod=xrtmethod, 
                         input_xrt_region=input_xrt_region, input_xrt_region_dict=input_xrt_region_dict,
                         xrt_exposure_dict=exposure_dict, use_prior_prep=use_prior_prep, aia_clobber=aia_clobber,
                            real_aia_err=real_aia_err, default_err=default_err, use_highTprep=True)

    #=========================================================================================================                 
    minT4=5.6
    maxT4=6.7
    
    if demmethod=='DEMREG':                  
        res4 = dodem(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, eis=False, dT=dT,
                    name=name, compare_fpm=False, combine_fpm=True,
                    nuenergies=nuenergies, nustar_path='./', pile_up_corr=pile_up_corr,
                    make_nustar=True, gtifile=gtifile, datapath=datapath, adjacent_grades=adjacent_grades,
                    regfile=regfile, edit_regfile=edit_regfile,  use_fit_regfile=use_fit_regfile,
                    COM_nustar_region=COM_nustar_region,  nuclobber=False, special_pha=special_pha,
                    xrt_path=xrt_path, plot=False, xrt_exclude=xrt_exclude, 
                    xrt_factor=xrt_factor, xrtmethod=xrtmethod, xrt_exposure_dict=exposure_dict,
                    input_xrt_region=input_xrt_region, input_xrt_region_dict=input_xrt_region_dict,
                    real_xrt_err=real_xrt_err,
                    aiamethod=aiamethod, input_aia_region=input_aia_region, aia_clobber=aia_clobber,
                    input_aia_region_dict=input_aia_region_dict,
                    just_prep=False,
                    gloci=1, mc_in=True, mc_rounds=100, 
                    plotresp=False, reg_tweak=reg_tweak, rgt_fact=rgt_fact, max_iter=max_iter, 
                    minT=minT4, maxT=maxT4, plotMK=plotMK, use_prior_prep=use_prior_prep,
                    real_aia_err=real_aia_err, default_err=default_err, use_highTprep=True) 
    
    
    if demmethod=='XRT_ITER':
        res4 = run_iterative_wrapper(time, bl, tr, minT4, maxT4, xrt=xrt, aia=aia, nustar=nustar, eis=False,
                          dT=dT, nuenergies=nuenergies, nustar_path='./', pile_up_corr=pile_up_corr,
                          edit_regfile=edit_regfile, use_fit_regfile=use_fit_regfile, adjacent_grades=adjacent_grades,
                          COM_nustar_region=COM_nustar_region, nuclobber=False, special_pha=special_pha,
                          xrt_factor=xrt_factor, real_xrt_err=real_xrt_err,
                          aia_exclude=aia_exclude, name=name,
                         plotMK=plotMK, mc_iter=100,
                         combine_fpm=True, aiamethod=aiamethod, input_aia_region=input_aia_region,
                         input_aia_region_dict=input_aia_region_dict, xrtmethod=xrtmethod, 
                         input_xrt_region=input_xrt_region, input_xrt_region_dict=input_xrt_region_dict,
                         xrt_exposure_dict=exposure_dict, use_prior_prep=use_prior_prep, aia_clobber=aia_clobber,
                            real_aia_err=real_aia_err, default_err=default_err, use_highTprep=True)
    
    #=========================================================================================================

    data1, timestring1 = vdr.load_DEM(time, 
                              filename=res1)
    data2, timestring2 = vdr.load_DEM(time,
                              filename=res2)
    data3, timestring3 = vdr.load_DEM(time, 
                              filename=res3)
    data4, timestring4 = vdr.load_DEM(time, 
                              filename=res4)
                       
    consistent = vdr.compare_DEMs(data1, data2, timestring1, timestring2, 
                                           title1=str(minT1)+'_'+str(maxT1)+'_DEM', 
                                           title2=str(minT2)+'_'+str(maxT2)+'_DEM')
    consistent = vdr.compare_DEMs(data1, data3, timestring1, timestring3,
                                           title1=str(minT1)+'_'+str(maxT1)+'_DEM', 
                                           title2=str(minT3)+'_'+str(maxT3)+'_DEM')
    consistent = vdr.compare_DEMs(data1, data4, timestring1, timestring4,
                                           title1=str(minT1)+'_'+str(maxT1)+'_DEM', 
                                           title2=str(minT4)+'_'+str(maxT4)+'_DEM')
                       
    
    return 
                       

def pickle_res(res, picklefile, plotMK, time, nuenergies, minT, maxT):
                       
    dem, dn_reg, ts, ts_, trmatrix, dn_in, edn_in, chanax, demerror, dnregerror, chisq, fpm = res
    data = {'time_interval': time,
            'nuenergies': nuenergies,
            'nufpm': fpm,
            'temp_interval': [minT, maxT],
            'plotMK': plotMK,
            'DEM': dem,
            'dn_reg': dn_reg,
            'ts': ts,
            'ts_': ts_,
            'trmatrix': trmatrix,
            'dn_in': dn_in,
            'edn_in': edn_in,
            'chanax': chanax,
            'edem': demerror,
            'edn': dnregerror, 
            'chisq': chisq
            }
    with open(picklefile, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)   
    
    
    return
    
    
    
    
    
    


