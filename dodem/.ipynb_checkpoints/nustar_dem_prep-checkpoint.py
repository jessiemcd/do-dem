import nustar_utilities as nuutil
import region_fitting as rf

import hissw

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.io as io

from os.path import splitext, isfile
import os 

from astropy.io import fits
from astropy import coordinates as coord
from astropy import units as u
import importlib
from nustar_pysolar import convert, utils
import nustar_pysolar.map as numap
from scipy import ndimage
import sunpy

import glob

import subprocess
import pathlib
import copy
import pickle

"""
Preparing NuSTAR Data for DEM: IDL/other helpers + other needed prep 
----------------------------------------------------------------------

        NuSTAR
        ------
        
        NuSTAR response to plasma is found by combining instrument response (from spectral data 
        products) and model photon spectrum from hot plasma as calculated via the fvth() model. 
        
        NuSTAR spectral data products for specific time/region:
        
                -If not found in timestring directory (named like /19-09-00_19-14-00/ after time 
                interval hhmmss_hhmmss), these are made.
                
        Model photon spectrum:
        
                -'fvth_out.dat' contains an array of fvth model photon spectrum.
                -If not present in working directory, this file is made by running make_fvtharray.pro 
                    through HISSW.
                    

        Steps to creating NuSTAR spectral data products (overview for humans – all functionality is 
        in make_nustar_products): 
        
            -Time Interval Selection:
                -Takes full .evt file from observation + input gti.fits file.
                -Makes new gti (good time interval) file for selected time interval.
                -Uses new gti file to create new .evt files with only data from that interval 
                (NUSCREEN pipeline step).
                -Does this for different event grades (0, 0-4, 21-24) to facilitate pileup correction (or not).
                
            -Region Validation:
                -Takes grade 0 .evt file and .reg file, converts .evt file to solar coordinates
                -Optimizes region (based on options chosen as inputs to rf.get_file_region)
                –Makes + plots NuSTAR map with region overlay, saves figure (to evaluate by eye if the 
                region is reasonable if needed at a later time).
                -Note that region selection is based on only the grade 0 events.
                
            –Generate Spectral Data Products:
                -Takes .evt files and .reg file
                -Makes spectral data products (NUPRODUCTS pipeline step). Does this for each grade. 
  

"""


def make_nustar_products(time, fpm, gtifile, datapath, regfile, nustar_path,
                         #region related inputs:
                         edit_regfile=True, centroid_region=False, nuradius=150,
                         twogauss=False, onegauss=False, direction='', guess=[], guess2=[],
                         #general method related:
                        compare_fpm=False, pile_up_corr=False, adjacent_grades=False,
                         clobber=False, path_to_dodem='/Users/jmdunca2/do-dem/', dip_before_products=False,
                        energy_percents=False):
    """
    See load_nustar() documentation.
    
    """
    
    
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
    #print(timestring)
    
    #Make time interval directory if it ddoessn't yet exist.
    save_path = pathlib.Path(nustar_path) / timestring
    if not save_path.exists():
        save_path.mkdir()

    grade_prep=0  
    if adjacent_grades:
        #USING GRADES 0-4
        arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, grade='0_4')
        if len(arf_files) == 1 and len(rmf_files) == 1 and len(pha_files) == 1:
            grade_prep+=1
        if pile_up_corr:
            arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, grade='21_24')
            if len(arf_files) == 1 and len(rmf_files) == 1 and len(pha_files) == 1:
                grade_prep+=1
            if grade_prep == 2:
                print('We have both grades 21-24 and grades 0-4 products already, and clobber is not set - exiting.')
                return
        else:
            if grade_prep == 1:
                print('We have grades 0-4 products already, and clobber is not set - exiting.')
                return
            
    else:
        #USING GRADE 0
        arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, grade='0')
        if len(arf_files) == 1 and len(rmf_files) == 1 and len(pha_files) == 1:
            grade_prep+=1
        
        if pile_up_corr:
            arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, grade='21_24')
            if len(arf_files) == 1 and len(rmf_files) == 1 and len(pha_files) == 1:
                grade_prep+=1
            if grade_prep == 2:
                print('We have both grades 21-24 and grade 0 products already, and clobber is not set - exiting.')
                return
        else:
            if grade_prep == 1:
                print('We have grade 0 products already, and clobber is not set - exiting.')
                return        
            
        
        
    # grade_prep=0    
    # arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, grade='0')
    # if len(arf_files) == 1 and len(rmf_files) == 1 and len(pha_files) == 1:
    #     grade_prep+=1
    # arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, grade='0_4')
    # if len(arf_files) == 1 and len(rmf_files) == 1 and len(pha_files) == 1:
    #     grade_prep+=1
    
    # if pile_up_corr==False:
    #     if grade_prep==2 and clobber==False:
    #         print('We have both grade 0 and grades 0-4 products already, and clobber is not set - exiting.')
    #         return
    # if pile_up_corr==True:
    #     arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, grade='21_24')
    #     if len(arf_files) == 1 and len(rmf_files) == 1 and len(pha_files) == 1:
    #         grade_prep+=1
    #     if grade_prep==3 and clobber==False:
    #         print('We have, grade 0, grades 0-4, and grades 21-24 products already, and clobber is not set - exiting.')
    #         return
    
    #======================================================
    #TIME INTERVAL SCREENING

    unphys_products=0
    if pile_up_corr:
        unphys_products=1

    if adjacent_grades:
        adj_value=1
    else:
        adj_value=0

    if adjacent_grades:
        evt_files_0 = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_4_p_cl.evt')
    else:
        evt_files_0 = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_p_cl.evt')

    #print(nustar_path+timestring+'/'+timestring+fpm+'_gti.fits')
    #print(nustar_path)
    #print(gtifile)
    
    if len(evt_files_0) != 1 or clobber==True:
        #If there is not an .evt files already for this FPM, (or clobber set), make a gti file for the 
        #time interval and fpm selected, and run nuscreen to make a time-interval-specific event list.
        edit_gti(gtifile, time[0], time[1], nustar_path+timestring+'/'+timestring+fpm+'_gti.fits')

        #print('prenuscreen')
        #print(path_to_dodem+'run_nuscreen.sh')
        #print(nustar_path+timestring+'/run_nuscreen.sh')
        #os.rename(path_to_dodem+'/run_nuscreen.sh', nustar_path+timestring+'/run_nuscreen.sh')
        status = subprocess.call('cp '+path_to_dodem+'run_nuscreen.sh '+nustar_path+timestring+'/run_nuscreen.sh', shell=True) 

        
        #Edit shell script to run nuscreen (component of NuSTAR pipeline). This makes grade 0, 0-4 and 21-24 .evt files.
        edit_nuscreen(nustar_path+timestring+'/', nustar_path, timestring, fpm, datapath, adjacent_grades=adj_value, 
                      unphys_products=unphys_products)
        f = open(nustar_path+timestring+'/'+fpm+"nuscreen_output.txt", "w")
        #screenprocess = subprocess.call(path_to_dodem+'run_nuscreen.sh', stdout=f)
        screenprocess = subprocess.run(nustar_path+timestring+'/run_nuscreen.sh', stdout=f, shell=True,
                                      cwd=nustar_path+timestring+'/')

        #print('postnuscreen')
        
    #Now we should definitely have the desired .evt file:    
    if adjacent_grades:
        evt_files_0 = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_4_p_cl.evt')
    else:
        evt_files_0 = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_p_cl.evt')

    #print('len evt files: ', len(evt_files_0))
        
    if len(evt_files_0) !=1:
        print('Failed to find or make grade 0 or grade 0-4 .evt files – not using NuSTAR.')
        return
    
    if pile_up_corr:
        evt_files_unphys = glob.glob(nustar_path+timestring+'/*'+fpm+'06_21_24_p_cl.evt')
        print('len evt files unphys: ', len(evt_files_unphys))
        if len(evt_files_unphys) !=1:
            print('Failed to find or make grade 21-24 .evt files, cant do pile-up correction. Not using NuSTAR.')
            return
    
    #======================================================

    #If we only wanted the time- and grade- specific evt files, get out now.

    if dip_before_products:
        print('dip_before_products is set True, so we will return after making only the time- and grade-specific evt files.')
        return

    #print('got here.')
    
    #======================================================
    #REGION SELECTION
 
    #Now, let's convert the grade 0 (or grade 0-4) .evt file to solar coordinates (if there isn't a solar coordinates 
    #file already:)
    
    if adjacent_grades:
        sun_file = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_4_p_cl_sunpos.evt')
        #If we don't already have the sunpos file...
        if len(sun_file) != 1 or clobber==True:
            evt_files_04 = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_4_p_cl.evt')
            try:
                convert_wrapper(evt_files_04[0], clobber=clobber)
            except (TimeoutError, OSError):
                print('Possible connection error (solar coordinate conversion requires internet) – process failed here.')
                return
            sun_file = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_4_p_cl_sunpos.evt')
            #If we still don't have the sunpos file:
            if len(sun_file) !=1:
                print('Failed to find or make grade 0-4 sunpos.evt file – not using NuSTAR.')
                return        
    else:
        sun_file = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_p_cl_sunpos.evt')
        #If we don't already have the sunpos file...
        if len(sun_file) != 1 or clobber==True:
            try:
                convert_wrapper(evt_files_0[0], clobber=clobber)
            except (TimeoutError, OSError):
                print('Possible connection error (solar coordinate conversion requires internet) – process failed here.')
                return
            sun_file = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_p_cl_sunpos.evt')
            #If we still don't have the sunpos file:
            if len(sun_file) !=1:
                print('Failed to find or make grade 0 sunpos.evt file – not using NuSTAR.')
                return

    #print('got past sunfiles.')

    #print('Twogauss set to: ', twogauss)
    if edit_regfile: 
        #print('Twogauss set to: ', twogauss)
        #Taking our solar-coordinates file, let's make a region in order to generate spectral data products!
        res = rf.get_file_region(sun_file[0], time[0], time[1], regfile, centroid_region=centroid_region, 
                                                 radius=nuradius, working_dir=nustar_path, efilter=False, 
                                                twogauss=twogauss, onegauss=onegauss, direction=direction, 
                                                 guess=guess, guess2=guess2)
        if res is None:
            return
        else:
            newregfile, percent = res

            
    else:
        newregfile=regfile
        plot_region=True
        if plot_region:
            rf.plot_file_region(sun_file[0], time[0], time[1], regfile)



    if energy_percents:
        percents(sun_file[0], newregfile, time, nustar_path+timestring+'/'+timestring+'_'+fpm+'_epercents.pickle')
        
    #======================================================
    
    #======================================================
    #SPECTRAL DATA PRODUCTS
    #Now that we have the region file + the time-interval-specific .evt file, let's do it!

    #print('prenuproducts')
    #os.popen('cp '+path_to_dodem+'/run_nuproducts.sh '+nustar_path+timestring+'/run_nuproducts.sh')
    status = subprocess.call('cp '+path_to_dodem+'/run_nuproducts.sh '+nustar_path+timestring+'/run_nuproducts.sh', shell=True) 
    
    #Edit shell script to run nuproducts (component of NuSTAR pipeline)
    edit_nuproducts(nustar_path+timestring+'/', nustar_path, timestring, fpm, newregfile, datapath, 
                    unphys_products=unphys_products,
                   adjacent_grades=adj_value)
    f = open(nustar_path+timestring+'/'+fpm+"nuproducts_output.txt", "w")
    productprocess = subprocess.run(nustar_path+timestring+'/run_nuproducts.sh', stdout=f, shell=True,
                                      cwd=nustar_path+timestring+'/')
    #print('postnuproducts')
    #======================================================

    if edit_regfile:
        return percent
    else:
        return newregfile
    
def find_nuproducts(nustar_path, timestring, fpm, special_pha='', grade='0', shush=False):
    """
    Looks at nustar_path + timestring directory for nustar spectral products for a given grade+fpm.
    Wrapper since we do this a few times.
    """
    arf_files = glob.glob(nustar_path+timestring+'/*'+fpm+'*'+grade+'_p_sr.arf')
    rmf_files = glob.glob(nustar_path+timestring+'/*'+fpm+'*'+grade+'_p_sr.rmf')
    #print(rmf_files)
    if bool(special_pha):
        pha_files = glob.glob(special_pha+timestring+'*'+fpm+'*.pha')
        if pha_files == []:
            print("Didn't find any .pha files in your special_pha directory.")
            print("Note expected format: timestring+'*'+fpm+'*.pha' where timestring is of the form 'hh-mm-ss_hh-mm-ss'")
    else:
        pha_files = glob.glob(nustar_path+timestring+'/*'+fpm+'06_'+grade+'_p_sr.pha')
    if not shush:
        print('ARF File: ', arf_files)
        print('RMF File: ', rmf_files)
        print('PHA File: ', pha_files)
    
    return arf_files, rmf_files, pha_files


def combine_fpm(time, eng_tr, nustar_path, make_nustar=False, gtifile='', datapath='', regfile='', 
                edit_regfile=True, actual_total_counts=False, centroid_region=False, use_fit_regfile=False,
               clobber=False, default_err=0.2, special_pha='', pile_up_corr=False, adjacent_grades=False,
               nuradius=150, path_to_dodem='./', countmin=10, force_both_fpm=False, shush=False,
                twogauss=False, onegauss=False, direction='', guess=[], guess2=[],
                   energy_percents=False):
    """
    LOADS BOTH FPM + ADDS TOGETHER THE RATES + RESPONSE. 
    
    See load_nustar() for documentation of inputs.
    
    """
    fpm = 'A'
    fpm2 = 'B'
    
    res = load_nustar(time, eng_tr, nustar_path, fpm, make_nustar=make_nustar,
                                          gtifile=gtifile, datapath=datapath, regfile=regfile, 
                                          edit_regfile=edit_regfile, compare_fpm=False,
                                             combine_fpm=True, actual_total_counts=actual_total_counts,
                                             centroid_region=centroid_region, use_fit_regfile=use_fit_regfile, clobber=clobber,
                                             default_err=default_err, special_pha=special_pha,
                                             pile_up_corr=pile_up_corr, adjacent_grades=adjacent_grades,
                                             nuradius=nuradius, path_to_dodem=path_to_dodem, shush=shush,
                                             twogauss=twogauss, onegauss=onegauss, direction=direction, guess=guess, guess2=guess2,
                                             energy_percents=energy_percents)
    if res is None:
        print('Something is wrong with ', fpm,'; Not using NuSTAR.')
        print('')
        if actual_total_counts:
            print('Region issue, returning 0 counts to force longer time interval in TIS.')
            return (0, False)
        else:  
            return
    
    if actual_total_counts:
        rate, erate, nutrs, nu_tresp, nu_logt, fpm, atc = res
        if atc >=countmin:
            print(atc, ' counts just in FPM', fpm)
            if force_both_fpm==True:
                print('Making products for FPM', fpm2, ' as well, as you set force_both_fpm==True.')
            else:
                print('Exiting.')
                return (atc, True)
    
    res2 = load_nustar(time, eng_tr, nustar_path, fpm2, make_nustar=make_nustar,
                                          gtifile=gtifile, datapath=datapath, regfile=regfile, 
                                          edit_regfile=edit_regfile, compare_fpm=False,
                                             combine_fpm=True, actual_total_counts=actual_total_counts,
                                             centroid_region=centroid_region, use_fit_regfile=use_fit_regfile, clobber=clobber,
                                              default_err=default_err, special_pha=special_pha,
                                             pile_up_corr=pile_up_corr, adjacent_grades=adjacent_grades,
                                             nuradius=nuradius, path_to_dodem=path_to_dodem, shush=shush,
                                             twogauss=twogauss, onegauss=onegauss, direction=direction, guess=guess, guess2=guess2,
                                              energy_percents=energy_percents)
    if res2 is None:
        print('Something is wrong with ', fpm,'; Not using NuSTAR.')
        print('')
        if actual_total_counts:
            print('Region issue, returning 0 counts to force longer time interval in TIS.')
            return (0, False)
        else:  
            return
    
    if actual_total_counts:
        rate, erate, nutrs, nu_tresp, nu_logt, fpm, atc = res
        rate2, erate2, nutrs2, nu_tresp2, nu_logt2, fpm2, atc2 = res2
        if atc+atc2 >= countmin:
            print(atc, ' counts in FPM', fpm, ' ', atc2, ' counts in FPM', fpm2,'. Total:', atc+atc2, ' Exiting.')
            return (atc+atc2, True)
        else:
            print('Not enough counts.')
            print('')
            return (atc+atc2, False)
            
    else:    
        rate, erate, nutrs, nu_tresp, nu_logt, fpm = res
        rate2, erate2, nutrs2, nu_tresp2, nu_logt2, fpm2 = res2
    
    
    #Add FPM to data labels
    fnutrs = [n+' '+fpm for n in nutrs]
    f2nutrs = [n+' '+fpm2 for n in nutrs2]
    #Number of energy ranges in FPM with the most ranges
    number = np.max([len(nutrs), len(nutrs2)])
    
    #If the first FPM has fewer energy ranges, add placeholders
    while len(fnutrs) < number:
        fnutrs.append('')
        rate = np.append(rate, 0)
        erate = np.append(erate, 0)
        nu_tresp = np.append(nu_tresp, np.zeros((len(nu_tresp),1)), axis=1)
       
    #If the second FPM has fewer energy ranges, add placeholders
    while len(f2nutrs) < number:
        f2nutrs.append('')
        rate2 = np.append(rate2, 0)
        erate2 = np.append(erate2, 0)
        nu_tresp2 = np.append(nu_tresp2, np.zeros((len(nu_tresp2),1)), axis=1)
    
    #For each energy range...
    trs=[]
    comboerate=[]
    for i in range(0, number):
        if i > len(nutrs2)-1:
            nutrs2.append('')
        if i > len(nutrs)-1:
            nutrs.append('')
        #print(i, len(nutrs2))
        #Make new combined label
        if nutrs[i] == nutrs2[i]:
            trs.append(nutrs[i]+' A+B')
        else:
            trs.append(fnutrs[i]+' + '+f2nutrs[i])
        
        comboerate.append(erate[i]+erate2[i])

    #New rate: simple sum
    comborate = rate + rate2  
    #New response: simple sum of responses
    combotresp = np.add(nu_tresp, nu_tresp2)
    
    
    return comborate, comboerate, trs, combotresp, nu_logt, fpm+'+'+fpm2
    
    


def load_nustar(time, eng_tr, nustar_path, fpm, make_nustar=False, gtifile='', datapath='', regfile='', 
                edit_regfile=True, compare_fpm=False, combine_fpm=False, actual_total_counts=False, centroid_region=False,
                   use_fit_regfile=False, clobber=False, default_err=0.2, pile_up_corr=False, special_pha='',
                   adjacent_grades=False, nuradius=150,path_to_dodem='./', shush=False,
                   twogauss=False, onegauss=False, direction='', guess=[], guess2=[], return_for_pile_up_figure=False,
                   energy_percents=False):
    """
    Load in NuSTAR data and response, return DEM inputs.
    
    Keywords
    ---------
    time - time interval for DEM (tuple of astropy Time objects)
                FORMAT LIKE, 
                time=(astropy.time.Time('2018-05-29T19:08:00', scale='utc'), 
                        astropy.time.Time('2018-05-29T19:14:00', scale='utc')
    
    eng_tr  - NuSTAR energy range(s) for DEM. 
                Format like eng_tr=[[2.5,3.5],[3.5,5], [5.,7.]]
                        or
                  eng_tr=[2.5,7.]
    
    nustar_path - location of NuSTAR spectral data products DIRECTORY (time interval directory) + 
                    fvth_out.dat (photon model). Or, where they will be placed if not found. 
    
    fpm - NuSTAR focal plane module (A or B). Used to identify spectral data products, etc.
    
    make_nustar - Set True to make NuSTAR spectral data products (if possible) for the 
                    indicated time interval, if they don't already exist in nustar_path.
                    
    clobber - Set True to re-make all the spectral data products (even if they already exist). 
    
    default_err - Factor defining what percent of the NuSTAR rates will be added in quadrature with counting
                    statistics to define the uncertainty on the DEM inputs (e.g. default_err=0.2 -> 20% 
                    uncertainty). 
    adjacent_grades - Set True to use grades 0-4 events as DEM inputs (if False, grade 0 will be used). 
               
    pile_up_corr - Set True to do a pileup correction. Current method: 
                    If adjacent_grades==True: (grades 0-4) - 5/4*(grades 21-24)
                    If adjacent_grades==False: (grade 0) - 1/4*(grades 21-24)
    
    special_pha - Set True to use pha files found at:  (special_pha+timestring+'*'+fpm+'*.pha') instead of
                    default format/those made by this code. Useful for making external modifications to 
                    spectrum (say, subtracting a non-thermal component), and then being able to plug a modified
                    pha file into this process.
                    
    path_to_dodem - path to location of do-dem installation (location of shell scripts for pipeline steps, as well 
                        as the hissw-idl directory).
                    
    Spectral Products Keywords
    --------------------------
                    
    gtifile - template file, will be changed based on time range
    
    datapath - path to NuSTAR data directory (top level, i.e. obsid-named directory)
    
    regfile - region file. Will be edited (optimal region found) if edit_regfile==True
    
    edit_regfile - See above
    
    centroid_region - set True to find optimal region via data center of mass (ignoring chip gap). If False, region will be found 
                via optimizing location to include the maximum amount of emission. 
    
    compare_fpm - Make both fpma, fpmb products, use whichever has more NuSTAR emission in its optimal region. 
    
    combine_fpm - when load_nustar is run from within combine_fpm, this is set True. This means that in a case where
                    there is no data in a given input energy range, a value of zero and a normal response will be 
                    returned for that energy range (rather than the range just not being included). 
                    
                    This allows the response to be accurate when combining FPMA, FPMB (just because there are no
                    counts in one FPM does not mean the response in that FPM should be excluded from the total 
                    response).
                    
    actual_total_counts - Set True to return the actual total number of counts in the highest NuSTAR energy range
                            (not livetime corrected, not a rate). This is added to be used when load_nustar() is
                            being used to optimize DEM time intervals (minimzing time interval duration while 
                            retaining > a certain # actual NuSTAR counts in each in the highest energy bin). 
                            
                            For more info, see find_intervals().
    
    """
    
    #If edit_regfile=True, this will be updated. 
    newregfile=[]
    

    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    timestring=timestring+'_'+stopstring
    
    if use_fit_regfile:
        #Overide
        edit_regfile=False
        
    
    if edit_regfile==False:
        if use_fit_regfile:
            print('Using previously-found fit region file:')
            regfile_ = glob.glob(nustar_path+timestring+'/*'+fpm+'06_cl_sunpos_fit_region.reg')
            regfile = regfile_[0]  
        else:
            print('Using input regfile:')
            print(regfile)
            
        
    #DEM interval duration
    dur = (time[1]-time[0]).to(u.s).value
    
    #Make time interval directory if it ddoessn't yet exist.
    save_path = pathlib.Path(nustar_path) / timestring
    if not save_path.exists():
        save_path.mkdir()
    
    #======================================================
    # Load in the pha, arf and rmf (NuSTAR RESPONSE)

    
    if adjacent_grades:
        print('Using grades 0-4 NuSTAR events.')
        arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm,
                                                        special_pha=special_pha, grade='0_4',
                                                         shush=shush)
    else:
        print('Using grade 0 NuSTAR events.')
        arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, 
                                                          special_pha=special_pha, grade='0',
                                                         shush=shush)
    
    if pile_up_corr:
        arf_files_unphys, rmf_files_unphys, pha_files_unphys = find_nuproducts(nustar_path, timestring, fpm,
                                                                                special_pha=special_pha, grade='21_24',
                                                                                 shush=shush)
       
    #Will be switched to True if any spectral data products are missing, or clobber is True.    
    remake_0=False
    remake_u=False
    
    #Try to load all:
    if len(arf_files) == 1 and len(rmf_files) == 1 and len(pha_files) == 1 and clobber==False:
        e_lo1, e_hi1, eff_area = nuutil.read_arf(arf_files[0])
        #print('NuSTAR effective area array (function of input photon energy, units cm^2), shape: ', eff_area.shape)
        e_lo2, e_hi2, rmf_mat = nuutil.read_rmf(rmf_files[0])
        #print('NuSTAR response matrix (input photon energy vs. PHA channel, units counts (unitless)), shape: ',
                                        #rmf_mat.shape)
        engs,cnts,lvtm,ontim=nuutil.read_pha(pha_files[0])
    else:
        remake_0=True
        
    if pile_up_corr:
        if len(arf_files_unphys) == 1 and len(rmf_files_unphys) == 1 and len(pha_files_unphys) == 1 and clobber==False:
            #loaded variables with _u added to name are DIFFERENT between grade 0 + unphysical grades, others identical.
            e_lo1, e_hi1, eff_area_u = nuutil.read_arf(arf_files_unphys[0])
            #print('NuSTAR effective area array (function of input photon energy, units cm^2), shape: ', eff_area.shape)
            e_lo2, e_hi2, rmf_mat_u = nuutil.read_rmf(rmf_files_unphys[0])
            #print('NuSTAR response matrix (input photon energy vs. PHA channel, units counts (unitless)), shape: ',
                                            #rmf_mat.shape)
            engs,cnts_u,lvtm,ontim=nuutil.read_pha(pha_files_unphys[0])
        else:
            remake_u=True
        
        

    #If we don't have all the response files...
    #Note: if pile_up_corr=True and there are not all the needed unphysical grade files, the grade 0 files will also be 
    #remade (and vice versa). However, if pile_up_corr=False, we will not make unphysical grade files.
    if remake_0 or remake_u:

        if not shush:
            print('There are ', len(arf_files), ' Auxiliary Response Files (.arf) for this FPM (', fpm, ') in nustar path.')
            print('There are ', len(rmf_files), ' Response Matrix Files (.rmf) for this FPM (', fpm, ') in nustar path.')
            print('There are ', len(pha_files), ' PHA files (.pha) for this FPM (', fpm, ') in nustar path.')
            print('We require 1 of each.')
            if pile_up_corr:
                print('There are ', len(arf_files_unphys), ' Unphysical Grade Auxiliary Response Files (.arf) for this FPM (', 
                      fpm, ') in nustar path.')
                print('There are ', len(rmf_files_unphys), ' Unphysical Grade Response Matrix Files (.rmf) for this FPM (', 
                      fpm, ') in nustar path.')
                print('There are ', len(pha_files_unphys), ' Unphysical Grade PHA files (.pha) for this FPM (', 
                      fpm, ') in nustar path.')
                print('We require 1 of each.')        
        
        if make_nustar:
            if compare_fpm:
                print('Comparing FPM, so we will edit the region file (setting edit_regfile=True)')
                edit_regfile=True
            if clobber:
                print('Clobber set– so we will make data products anyway.')
            else:
                print('Now we will make some spectral data products.')
            mn = make_nustar_products(time, fpm, gtifile, datapath, regfile, nustar_path, edit_regfile=edit_regfile, 
                                      centroid_region=centroid_region, clobber=True, pile_up_corr=pile_up_corr, 
                                      adjacent_grades=adjacent_grades, nuradius=nuradius,
                                     twogauss=twogauss, onegauss=onegauss, direction=direction, guess=guess, guess2=guess2,
                                     energy_percents=energy_percents)
            if compare_fpm:
                if fpm == 'A':
                    fpm2='B'
                else:
                    fpm2='A'
                print('Compare FPM is set – examining fpm', fpm2, 
                          ' also to see which has more emission in its optimal region.')
                mn2 = make_nustar_products(time, fpm2, gtifile, datapath, regfile, nustar_path, edit_regfile=edit_regfile, 
                                      centroid_region=centroid_region, clobber=True, pile_up_corr=pile_up_corr, 
                                           adjacent_grades=adjacent_grades, nuradius=nuradius,
                                         twogauss=twogauss, onegauss=onegauss, direction=direction, guess=guess, guess2=guess2,
                                          energy_percents=energy_percents)
                print('FPM', fpm, ' region has ', mn, '% of emission. FPM', fpm2, ' region has ', mn2, '% of emission.' )
                if mn2>mn:
                    print('Switching to FPM', fpm2)
                    fpm=fpm2
            else:
                newregfile=mn
                
            if mn is None:
                return
            else:
                if adjacent_grades:
                    arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, grade='0_4')
                else:
                    arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, grade='0')

                if len(arf_files) == 1 and len(rmf_files) == 1 and len(pha_files) == 1:
                    print('now we have everything for physical grades.')
                else:
                    print('still failed to make all physical grade files.')
                    return
                    
                e_lo1, e_hi1, eff_area = nuutil.read_arf(arf_files[0])
                e_lo2, e_hi2, rmf_mat = nuutil.read_rmf(rmf_files[0])
                engs,cnts,lvtm,ontim=nuutil.read_pha(pha_files[0])
                if pile_up_corr:
                    arf_files_unphys, rmf_files_unphys, pha_files_unphys = find_nuproducts(nustar_path, timestring, fpm,
                                                                               special_pha=special_pha, grade='21_24')

                    if len(arf_files_unphys) == 1 and len(rmf_files_unphys) == 1 and len(pha_files_unphys) == 1:
                        print('now we have everything for unphysical grades.')
                    else:
                        print('still failed to make all unphysical grade files.')
                        return
                    
                    e_lo1, e_hi1, eff_area_u = nuutil.read_arf(arf_files_unphys[0])
                    e_lo2, e_hi2, rmf_mat_u = nuutil.read_rmf(rmf_files_unphys[0])
                    engs,cnts_u,lvtm,ontim=nuutil.read_pha(pha_files_unphys[0]) 
                    
                    
        else:
            print('You need to set make_nustar=True to make spectral data products.')
            return
    
    print('')
    
    #Get area of NuSTAR region
    if centroid_region or twogauss or onegauss:
        area = nuradius**2*np.pi*7.25e7*7.25e7  
    else:
        if bool(newregfile)==False:
            newregfile=regfile
        print("For NuSTAR area, using region in: ", newregfile)
        offset, rad = rf.read_regfile(newregfile, time[0], time[1], 'hourangle')
        area = (rad.value)**2*np.pi*7.25e7*7.25e7 


    #======================================================
    
    
    #======================================================
    #LOAD IN THE THERMAL EMISSION MODEL 

    #Output of make_fvtharray.pro (sswidl)
    try:
        fvth=io.readsav('fvth_out.dat')
    except FileNotFoundError:
        print('Hey, where is the thermal emission model file? It should be called:')
        print(nustar_path+'fvth_out.dat')
        print('Let us now make a new one using hissw to run the IDL code.')
        
        ssw = hissw.Environment(ssw_packages=['xray', 'hessi'])
        agr_path = path_to_dodem+'/hissw_idl/fvtharray_hissw_wrapper.pro'
        try:
            ssw_resp = ssw.run(agr_path)
            fvth=io.readsav('fvth_out.dat')
        except Exception:
            import traceback
            print(traceback.print_exc())
            print('Something is wrong with the SSWIDL run - make sure the following IDL script exists:')
            print(agr_path)
            print('')
            return
        
    #======================================================
    
    #======================================================
    #PREP DATA + RESPONSE FOR EACH ENERGY RANGE

    #List of energy bins
    engs_=fvth['eng'] #units keV
    #energy bin size
    de=engs_[1]-engs_[0]    
    #Only the part of the spectrum for which we have a thermal emission model
    cut_cnts=cnts[0:len(engs_)]
    
    
    #indices where there are counts
    wnz = np.where(cut_cnts > 0.)[0]
    enz = engs_[wnz]
    #ccnz = cut_cnts[wnz]
    
    if pile_up_corr:
        fig=plt.figure(figsize=(10,5))
        
        if adjacent_grades:
            plt.plot(engs, (cnts-(5./4)*cnts_u), label='corr')
            plt.plot(engs, cnts, label='og grades 0-4')
            plt.plot(engs, cnts_u, label='unphys')
            plt.yscale('log')
            plt.xlim([0,20])
            plt.legend()
            plt.savefig(nustar_path+timestring+'/'+timestring+fpm+'_adjacent_pile_up.png')   
            plt.close(fig)
            
            cnts_corr = cnts-(5./4)*cnts_u
        else:
            plt.plot(engs, (cnts-0.25*cnts_u), label='corr')
            plt.plot(engs, cnts, label='og grade 0')
            plt.plot(engs, cnts_u, label='unphys')
            plt.yscale('log')
            plt.xlim([0,20])
            plt.legend()
            plt.savefig(nustar_path+timestring+'/'+timestring+fpm+'pile_up.png')
            plt.close(fig)
            
            cnts_corr = cnts-0.25*cnts_u

        if return_for_pile_up_figure:
            return engs, cnts, cnts_u
    
    numax = enz[-1]
    
    #print('Max NuSTAR Energy: ', numax)
    if numax > eng_tr[-1][1]:
        print('Warning: there is at least one NuSTAR event of higher energy than your highest energy range')
        print('Highest energy range:', eng_tr[-1])
        
        highEs = enz[np.where(enz > eng_tr[-1][1])[0]]
        print('Total Above: ', len(highEs))
        print('Above Energies: ', highEs)
        
        #Uncomment to switch to changing the highest energy range to include all emission.
        #newtop = math.ceil(numax)
        #eng_tr[-1] = [eng_tr[-1][0], float(newtop)]
        #print('Changing highest energy range to:', eng_tr[-1])
    
    if combine_fpm==False:
        while numax < eng_tr[-1][0]:
            print('Warning: there are no counts in the highest energy range, ', eng_tr[-1])
            if len(eng_tr) == 1:
                print('AKA no counts above our minimum energy of interest– not using NuSTAR.')
                return
            eng_tr=eng_tr[0:-1]
            print('Switching to use only ', eng_tr)

    #List of temperatures (log space)
    logt=fvth['logt'] #units log(K)

    #Photon model: for each temperature, an array of what photon energies would be emitted.
        #Details: at each temperature, we find the array of counts emitted at each photon energy, 
        #assuming a isothermal spectral model with EM=1e49 cm^-3 and coronal abundances. Energy 
        #bin size: 0.04 keV.
    phmod=np.array(fvth['fvth']) # in units of photons/s/keV/cm2
    #Set nan (null) values to be 0, to keep them from messing things up later.
    phmod[np.isnan(phmod)] = 0.
    

    #Number of energies in f_vth output:
    nume=len(engs_)
    #The effective energy bins from the NuSTAR ARF, but only up to the highest energy bin from fvth.
    #i.e. effective area at each energy observed
    arf=eff_area[:nume] #units cm^2
    #segment of rmf matrix from the NuSTAR RMF, from 0 to the highest energy bin from fvth in each dimension
    rmf=rmf_mat[0:nume,0:nume] #units counts/photon


    #SRM Matrix:
    #Multiply the effective area at each input photon energy (value) by the response of the system to that 
    #input photon energy (array of PHA channels), to get a matrix relating input photon energy to measured
    #PHA channel that includes effective area information. 

    srm = np.array([rmf[r, :] * arf[r] for r in range(len(arf))]) #units cm^2*counts/photon
    
    #Uncomment to extract the response + spectrum and pickle it to examine elsewhere.
#     data = {'srm': srm,
#            'engs': engs,
#            'cnts': cnts}

#     with open('heres_an_srm.pickle', 'wb') as f:
#         # Pickle the 'data' dictionary using the highest protocol available.
#         pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)        

    # Defining two values corresponding to the dimensions of the photon model: n1 is the # of photon energies 
    #in the photon model, n2 is the # of temperatures.
    n1,n2=phmod.shape
    # Making a zero array with the shape of the photon model 
    modrs= np.zeros([n1,n2])

    # For each temperature in the photon model, make a row in this new matrix consisting of:
    # the array of counts in each photon energy bin generated by isothermal plasma at that temperature
    # multiplied BY
    # the SRM matrix (input photon vs. recorded PHA channel)
    # giving us:
    # a new "Model Response Matrix" relating the plasma temperature to the observed NuSTAR PHA counts
    for t in np.arange(n2):
        modrs[:,t]=(phmod[:,t]@srm)*de #photons/s/keV/cm2 * cm^2*counts/photon * keV = counts/s

    #Make zeros array with dimensions of: number of temperatures, number of energy ranges
    tresp=np.zeros([len(modrs[0,:]),len(eng_tr)])
    #print(tresp.shape)

    for i in np.arange(len(eng_tr)):
        erange=eng_tr[i]
        #indices where e_lo is greater than the low end of the energy range, and e_hi is less than the high end
        #AKA energy bins within the energy range
        gd=np.where((e_lo1 >= erange[0]) & (e_hi1 < erange[1]) )

        #mm: vector with sum of the counts in all those energy bins at each temperature
        mm=np.sum(modrs[gd,:], axis=1) #units counts/s

        #temperature response vector: expected NuSTAR counts in given energy range as a function
        #of temperature, divided by a factor of 10^49 (the input emission measure)

        #Emission measure acts as a scaling factor on spectrum: counts as function of energy increase/decrease but
        #spectral shape remains the same. Thus, dividing out the input emission measure makes this a kind of 
        #EM-independent observed count rate as a function of temp.
        
        #The input to fvth is unitless (IDL) - it's just 1 (input should be "em_49, emission measure units of 10^49"
        #according to documentation). Assuming volume EM units of 1/cm^3 here. Then, we multiply the response by the 
        #area used for data selection to give overall units counts/s * cm^5.

        tresp[:,i]=mm[0,:]*area/1e49 #units counts/s * cm^2 * cm^3

        
    # Work out the total count rate and error in the energy bands
    rate=np.zeros(len(eng_tr))
    if pile_up_corr:
        rate_corr=np.zeros(len(eng_tr))
    erate=np.zeros(len(eng_tr))
    if pile_up_corr:
        erate_corr=np.zeros(len(eng_tr))
    nutrs=[]
        
    #print('E Rates NuSTAR')

    quitindex=[]
    quit=False
    for i in np.arange(len(eng_tr)):
        erange=eng_tr[i]
        #Where the NuSTAR-observed energies (from PHA file) were within the selected energy range.
        gd=np.where((engs >= erange[0]) & (engs < erange[1]) )
        #Rate: (units:cts/s) energies in energy range divided by the livetime
        if pile_up_corr:
            ratt = np.sum(cnts[gd])/lvtm
            rate[i]=ratt
            ratt_c = np.sum(cnts_corr[gd])/lvtm
            if ratt_c < 0:
                quitindex=i
                quit=True
                continue
            rate_corr[i]=ratt_c
            erate_corr[i]=np.sqrt(np.sum(cnts_corr[gd]))/lvtm
            #print('Erate, corr:', erate_corr[i])
            #print('Default % of rate, squared:', (default_err*rate_corr[i])**2)
            erate_corr[i]=(erate_corr[i]**2+(default_err*rate_corr[i])**2)**0.5
            
            if actual_total_counts:
                atc = np.sum(cnts_corr[gd])
                print('ATC:', atc)
                print('')
            
            
        else:
            ratt = np.sum(cnts[gd])/lvtm
            rate[i]=ratt
            #Error on the rate: square root of the counts divided by the livetime
            erate[i]=np.sqrt(np.sum(cnts[gd]))/lvtm
            #add previous error in quadrature with 20% uncertainty
            erate[i]=(erate[i]**2+(default_err*rate[i])**2)**0.5
            
            if actual_total_counts:
                atc = np.sum(cnts[gd])
                print('ATC:', atc)
                print('')
        
        #Make label
        nutr=str(erange[0])+'-'+str(erange[1])+'keV'
        nutrs.append(nutr)
       
            
    #======================================================
    
    if pile_up_corr:
        #Uncomment to check on magnitude of pileup correction
        #for r in range(0,len(rate_corr)):
            #print('Pileup Correction in ', eng_tr[r], ': ', round((rate[r]-rate_corr[r])/rate[r]*100,3), '%')
            #if rate[r] < rate_corr[r]:
            #    print('Pileup Correction Added Counts?!: og - ', rate[r], ' corr - ', rate_corr[r])
        #print('Rates:', rate)
        #print('Corrected Rates:', rate_corr)
        #print('Corrected ERates:', erate_corr)
            
        rate=rate_corr
        erate=erate_corr
       
    if bool(quit):
        print('')
        print(nutrs)
        print('Removed rate in ', eng_tr[quitindex], fpm, ' for being negative after pile-up correction.')
        rate = rate[0:quitindex]
        erate = erate[0:quitindex]
        tresp = tresp[:, 0:quitindex]   
        atc = 0
    
    if actual_total_counts:
        return rate, erate, nutrs, tresp, logt, fpm, atc
    else:
        return rate, erate, nutrs, tresp, logt, fpm


def percents(evt_file, regfile, time_interval, savefile):

    import nustar_pysolar as nustar
    from regions import CircleSkyRegion

    #print(evt_file)
    
    with fits.open(evt_file) as hdu:
        evt_data = hdu[1].data
        hdr = hdu[1].header
    
    bounds=[[0.,20.],[2.5,20.],[2.5, 3.5], [3.5,6.], [6.,10.], [10.,20.]]#,10.]

    percentdict = {'evt_file': evt_file,
                  'reg_file': regfile}
    
    
    for b in bounds:
        perstring = str(b[0])+'-'+str(b[1])+' keV'
        
        cleanevt = nustar.filter.event_filter(evt_data, energy_low=b[0], energy_high=b[1],
                                             no_bad_pix_filter=True, no_grade_filter=True)
        #print(len(cleanevt), len(evt_data))
        try:
            nustar_map = nustar.map.make_sunpy(cleanevt, hdr)
        except ValueError:
            print('Failed on range ', b, ' keV - no info')
            percentdict[perstring] = float("nan")
            continue
    
    
        offset, rad = rf.read_regfile(regfile, time_interval[0], time_interval[1], 'hourangle')
        region = CircleSkyRegion(
                        center = coord.SkyCoord(offset[0], offset[1], frame=nustar_map.coordinate_frame),
                        radius = rad
                    )
        
        regdata = rf.get_region_data(nustar_map, region, 0)
        percent = np.sum(regdata)/np.sum(nustar_map.data)
        print('Percent of emission between '+str(b[0])+', '+str(b[1])+' keV in region:', round(percent*100,1))

        percentdict[perstring] = percent


    with open(savefile, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(percentdict, f, pickle.HIGHEST_PROTOCOL) 



def edit_nuscreen(path_to_run_nuscreen, nustar_path, timestring, fpm, datapath, adjacent_grades=0, 
                      unphys_products=0):
    """
    Requires shell script run_nuscreen.sh in input path (nustar_path). 
    Edits based on inputs + overwrites new version.
    """

    with open(path_to_run_nuscreen+'run_nuscreen.sh', "r+") as f:
        screen = f.read()
        #print(screen)

        screenlist = screen.split('\n')
        #print(screenlist)
        screenlist[4] = 'interval='+timestring
        screenlist[5] = 'fpm='+fpm
        screenlist[6] = 'INDIR='+datapath
        screenlist[7] = 'work_dir='+nustar_path
        screenlist[8] = 'unphys_products='+str(unphys_products)
        screenlist[9] = 'adjacent_grades='+str(adjacent_grades)

        newscreen = '\n'.join(screenlist)
        f.seek(0)
        f.write(newscreen)
        f.truncate()
            
            
def edit_nuproducts(path_to_run_nuproducts, nustar_path, timestring, fpm, regfile, datapath, unphys_products=0,
                   adjacent_grades=0):
    """
    Requires shell script run_nuproducts.sh in input path (nustar_path). 
    Edits based on inputs + overwrites new version.
    """
    
    with open(path_to_run_nuproducts+'run_nuproducts.sh', "r+") as f:
        nuproducts = f.read()
        productslist = nuproducts.split('\n')
        #print(productslist)
        productslist[4] = 'interval='+timestring
        #productslist[5] = 'region='+regfile.split('/')[-1][0:-4]
        productslist[5] = 'region='+regfile
        productslist[6] = 'fpm='+fpm
        productslist[7] = 'INDIR='+datapath
        productslist[8] = 'unphys_products='+str(unphys_products)
        productslist[9] = 'working_dir='+nustar_path
        productslist[10] = 'adjacent_grades='+str(adjacent_grades)

        newnuproducts = '\n'.join(productslist)
        #print(newnuproducts)
        f.seek(0)
        f.write(newnuproducts)
        f.truncate()


def convert_wrapper(infile, clobber=False):
    """
    Personalized wrapper for using nustar_pysolar to convert an .evt file (infile) to solar coordinates.
    
    (Needed to make images to check region file reasonability before making spectral data products).
    
    """
    
    # # Make the new filename:
    (sfile, ext)=splitext(infile)
    outfile=sfile+'_sunpos.evt'

    if clobber==False:
        # If we already have one of the sunpos files, don't re-do that one too.
        if isfile(outfile):
            #print(outfile, 'exists! Not re-making it.')
            return

    hdulist = fits.open(infile)
    evtdata = hdulist[1].data 
    hdr = hdulist[1].header
    hdulist.close()

    importlib.reload(convert)
    (newdata, newhdr) = convert.to_solar(evtdata, hdr)
    fits.writeto(outfile, newdata, newhdr, overwrite=clobber)
    
    return
    



def edit_gti(gti_file, newstart, newstop, newfile):
    """
    Takes existing NuSTAR GTI file and edits the start and stop times for the "good time interval".
    Please take care to use an informed new time interval (e.g. actually during the orbit of interest.)
    
    If a gti file with the same name as newfile already exists, it will be overwritten.

    Keywords
    ---------

    gti_file - path to old gti_file
    
    newfile - path/name for new gti file

    newstart, newstop = astropy.time.Time objects

    """
    from astropy.io import fits

    nunewstart = nuutil.convert_nustar_time(newstart, from_datetime=True, astropy_time=True)
    nunewend = nuutil.convert_nustar_time(newstop, from_datetime=True, astropy_time=True)
    with fits.open(gti_file) as hdu:
        hdr = hdu[1].header
        d = hdu[1].data
        d['START'] = np.array([round(nunewstart)])
        d['STOP'] = np.array([round(nunewend)])
        #print(hdu[1].data)
        
    fits.writeto(newfile, d, header=hdr, overwrite=True)


    

def plot_grade_spectra(working_dir, timestring, fpm):

    """
    Wrapper to plot spectra.
    """
    
    #unphysical grades
    arf_files_unphys, rmf_files_unphys, pha_files_unphys = find_nuproducts(working_dir, timestring, fpm,
                                                                               special_pha=False, grade='21_24')
    engs,cnts_u,lvtm,ontim=nuutil.read_pha(pha_files_unphys[0]) 

    #Grade 0 files: 

    arf_files, rmf_files, pha_files = find_nuproducts(working_dir, timestring, fpm, grade='0')
    engs,cnts,lvtm,ontim=nuutil.read_pha(pha_files[0])
    
    fig=plt.figure(figsize=(10,5))
    plt.plot(engs, (cnts-0.25*cnts_u), label='corr')
    plt.plot(engs, cnts, label='og grade 0')
    plt.plot(engs, cnts_u, label='unphys')
    plt.yscale('log')
    plt.xlim([0,13])
    plt.legend()
    plt.savefig(working_dir+timestring+'/'+timestring+fpm+'pile_up.png')
    plt.close(fig)


    #Grade 0-4 files: 
    
    arf_files, rmf_files, pha_files = find_nuproducts(working_dir, timestring, fpm, grade='0')
    engs,cnts,lvtm,ontim=nuutil.read_pha(pha_files[0])
    
    fig=plt.figure(figsize=(10,5))
    plt.plot(engs, (cnts-(5./4)*cnts_u), label='corr')
    plt.plot(engs, cnts, label='og grades 0-4')
    plt.plot(engs, cnts_u, label='unphys')
    plt.yscale('log')
    plt.xlim([0,13])
    plt.legend()
    plt.savefig(working_dir+timestring+'/'+timestring+fpm+'_adjacent_pile_up.png')   
    plt.close(fig)



def return_submap(datapath='./', fpm='A', specific_evt=[], 
                  bl=[], tr=[], nusmooth=True,
                  return_evt_hdr=False, plot=False,
                     return_obsid=False, already_sunpos=False):
    """
    wrapper - convert to solar coordinates and make submap for nice plot
    """
    from astropy.coordinates import SkyCoord
    
    #Get evt file
    if specific_evt:
        evt_file=specific_evt
    else:
        #print(datapath)
        evt_file = glob.glob(datapath+'/event_cl/*'+fpm+'06_cl.evt')[0]

    if already_sunpos:
        sun_file = evt_file
    else:
        #print(evt_file)
        #Convert evt file to solar coordinates         
        convert_wrapper(evt_file, clobber=False)
    
        #Get solar coordinates file
        if specific_evt:
            sun_file = specific_evt[:-4]+'_sunpos.evt'
        else:
            sun_file = glob.glob(datapath+'/event_cl/*'+fpm+'06_cl_sunpos.evt')[0]
        #print('Using solar coodinates file:', sun_file)
    
    with fits.open(sun_file) as hdu:
        evt_data = hdu[1].data
        hdr = hdu[1].header

    

    if return_evt_hdr:
        if return_obsid:
            obsid = sun_file.split('/')[-1][2:13]
            return evt_data, hdr, obsid
        return evt_data, hdr


    nustar_map = numap.make_sunpy(evt_data, hdr, norm_map=True)
    
    if nusmooth:
       from scipy import ndimage
       import sunpy.map
       #Smoothing the data; change sigma to smooth more or less
       dd=ndimage.gaussian_filter(nustar_map.data, sigma=2, mode='nearest')
       nustar_map=sunpy.map.Map(dd, nustar_map.meta)
        
    if bl:
        submap = nustar_map.submap(bottom_left=SkyCoord(*bl, frame=nustar_map.coordinate_frame),
                          top_right=SkyCoord(*tr, frame=nustar_map.coordinate_frame))
    else:    
        bl = SkyCoord( *(-1250, -1250)*u.arcsec, frame=nustar_map.coordinate_frame)
        tr = SkyCoord( *(1250, 1250)*u.arcsec, frame=nustar_map.coordinate_frame)
        submap = nustar_map.submap(bottom_left=bl, top_right=tr)

    if plot:
        print(submap.unit)
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(1,1,1, projection=submap)
        submap.plot(axes=ax)
        submap.draw_contours(5, axes=ax)#, index=1, percent=True, fill=True)

        #print(submap.contour(0.05))
        #ax.plot_coord(submap.contour(0.05))
        #print(submap.contour(0.05))

    #print(np.max(submap.data))
    #print(np.min(submap.data))
    #print(np.mean(submap.data))

    return submap
