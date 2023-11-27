import nustar_utilities as nuutil
import region_fitting as rf

import hissw

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.io as io

from os.path import splitext, isfile
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


def make_nustar_products(time, fpm, gtifile, datapath, regfile, nustar_path, edit_regfile=True,
                        compare_fpm=False, nofit=False, pile_up_corr=False, clobber=False, nuradius=150,
                        path_to_dodem='./'):
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
    
    #======================================================
    #TIME INTERVAL SCREENING
    
    evt_files_0 = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_p_cl.evt')
    
    if len(evt_files_0) != 1 or clobber==True:
        #If there is not an .evt files already for this FPM, (or clobber set), make a gti file for the 
        #time interval and fpm selected, and run nuscreen to make a time-interval-specific event list.
        edit_gti(gtifile, time[0], time[1], nustar_path+timestring+'/'+timestring+fpm+'_gti.fits')

        #Edit shell script to run nuscreen (component of NuSTAR pipeline). This makes both grade 0 and 21-24 .evt files.
        edit_nuscreen(path_to_dodem, nustar_path, timestring, fpm, datapath)
        f = open(nustar_path+timestring+'/'+fpm+"nuscreen_output.txt", "w")
        screenprocess = subprocess.call(path_to_dodem+'run_nuscreen.sh', stdout=f)
        
    #Now we should definitely have the desired .evt file:    
    evt_files_0 = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_p_cl.evt')
    if len(evt_files_0) !=1:
        print('Failed to find or make grade 0 .evt files – not using NuSTAR.')
        return
    
    if pile_up_corr:
        evt_files_unphys = glob.glob(nustar_path+timestring+'/*'+fpm+'06_21_24_p_cl.evt')
        if len(evt_files_unphys) !=1:
            print('Failed to find or make grade 21-24 .evt files, cant do pile-up correction. Not using NuSTAR.')
            return
    
    #======================================================
    
    #======================================================
    #REGION SELECTION
 
    #Now, let's convert the grade 0 .evt file to solar coordinates (if there isn't a solar coordinates 
    #file already:)
    sun_file_0 = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_p_cl_sunpos.evt')
    #If we don't already have the sunpos file...
    if len(sun_file_0) != 1 or clobber==True:
        convert_wrapper(evt_files_0[0], clobber=clobber)
        sun_file_0 = glob.glob(nustar_path+timestring+'/*'+fpm+'06_0_p_cl_sunpos.evt')
        #If we still don't have the sunpos file:
        if len(sun_file_0) !=1:
            print('Failed to find or make grade 0 sunpos.evt file – not using NuSTAR.')
            return
        
    if edit_regfile: 
        #Taking our solar-coordinates file, let's make a region in order to generate spectral data products!
        newregfile, percent = rf.get_file_region(sun_file_0[0], time[0], time[1], regfile, nofit=nofit, 
                                                 radius=nuradius,working_dir=nustar_path)
    else:
        newregfile=regfile
        
    #======================================================
    
    #======================================================
    #SPECTRAL DATA PRODUCTS
    #Now that we have the region file + the time-interval-specific .evt file, let's do it!
    unphys_products=0
    if pile_up_corr:
        unphys_products=1
        
    #Edit shell script to run nuproducts (component of NuSTAR pipeline)
    edit_nuproducts(path_to_dodem, nustar_path, timestring, fpm, newregfile, datapath, unphys_products=unphys_products)
    f = open(nustar_path+timestring+'/'+fpm+"nuproducts_output.txt", "w")
    productprocess = subprocess.call(path_to_dodem+'run_nuproducts.sh', stdout=f)
    
    #======================================================

    if edit_regfile:
        return percent
    else:
        return newregfile
    
def find_nuproducts(nustar_path, timestring, fpm, special_pha='', grade='0'):
    """
    Looks at nustar_path + timestring directory for nustar spectral products for a given grade+fpm.
    Wrapper since we do this a few times.
    """
    arf_files = glob.glob(nustar_path+timestring+'/*'+fpm+'*'+grade+'_p_sr.arf')
    rmf_files = glob.glob(nustar_path+timestring+'/*'+fpm+'*'+grade+'_p_sr.rmf')
    if bool(special_pha):
        pha_files = glob.glob(special_pha+timestring+'*'+fpm+'*.pha')
        if pha_files == []:
            print("Didn't find any .pha files in your special_pha directory.")
            print("Note expected format: timestring+'*'+fpm+'*.pha' where timestring is of the form 'hh-mm-ss_hh-mm-ss'")
    else:
        pha_files = glob.glob(nustar_path+timestring+'/*'+fpm+'06_'+grade+'_p_sr.pha')
    print('ARF File: ', arf_files)
    print('RMF File: ', rmf_files)
    print('PHA File: ', pha_files)
    
    return arf_files, rmf_files, pha_files


def combine_fpm(time, eng_tr, nustar_path, make_nustar=False, gtifile='', datapath='', regfile='', 
                edit_regfile=True, actual_total_counts=False, nofit=False, use_fit_regfile=False,
               clobber=False, default_err=0.2, special_pha='', pile_up_corr=False, adjacent_grades=False,
               nuradius=150, path_to_dodem='./'):
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
                                             nofit=nofit, use_fit_regfile=use_fit_regfile, clobber=clobber,
                                             default_err=default_err, special_pha=special_pha,
                                             pile_up_corr=pile_up_corr, adjacent_grades=adjacent_grades,
                                             nuradius=nuradius, path_to_dodem=path_to_dodem)
    if res is None:
        print('Something is wrong with ', fpm,'; Not using NuSTAR.')
        print('')
        return
    
    if actual_total_counts:
        rate, erate, nutrs, nu_tresp, nu_logt, fpm, atc = res
        if atc >=10:
            print(atc, ' counts just in FPM', fpm, '. Exiting.')
            return (atc, True)
    
    res2 = load_nustar(time, eng_tr, nustar_path, fpm2, make_nustar=make_nustar,
                                          gtifile=gtifile, datapath=datapath, regfile=regfile, 
                                          edit_regfile=edit_regfile, compare_fpm=False,
                                             combine_fpm=True, actual_total_counts=actual_total_counts,
                                             nofit=nofit, use_fit_regfile=use_fit_regfile, clobber=clobber,
                                              default_err=default_err, special_pha=special_pha,
                                             pile_up_corr=pile_up_corr, adjacent_grades=adjacent_grades,
                                             nuradius=nuradius, path_to_dodem=path_to_dodem)
    if res2 is None:
        print('Something is wrong with ', fpm,'; Not using NuSTAR.')
        print('')
        return
    
    if actual_total_counts:
        rate, erate, nutrs, nu_tresp, nu_logt, fpm, atc = res
        rate2, erate2, nutrs2, nu_tresp2, nu_logt2, fpm2, atc2 = res2
        if atc+atc2 >= 10:
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
                edit_regfile=True, compare_fpm=False, combine_fpm=False, actual_total_counts=False, nofit=False,
                   use_fit_regfile=False, clobber=False, default_err=0.2, pile_up_corr=False, special_pha='',
                   adjacent_grades=False, nuradius=150,path_to_dodem='./'):
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
    
    nofit - set True to find optimal region via data center of mass (ignoring chip gap). If False, region will be found 
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
                            retaining >10 actual NuSTAR counts in each in the highest energy bin). 
                            
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
    
    arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm, special_pha=special_pha, grade='0')
    
    if adjacent_grades:
            print('Using grades 0-4 NuSTAR events.')
            arf_files, rmf_files, pha_files = find_nuproducts(nustar_path, timestring, fpm,
                                                                               special_pha=special_pha, grade='0_4')
    
    if pile_up_corr:
        arf_files_unphys, rmf_files_unphys, pha_files_unphys = find_nuproducts(nustar_path, timestring, fpm,
                                                                               special_pha=special_pha, grade='21_24')
       
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
                                      nofit=nofit, clobber=True, pile_up_corr=pile_up_corr, nuradius=nuradius)
            if compare_fpm:
                if fpm == 'A':
                    fpm2='B'
                else:
                    fpm2='A'
                print('Compare FPM is set – examining fpm', fpm2, 
                          ' also to see which has more emission in its optimal region.')
                mn2 = make_nustar_products(time, fpm2, gtifile, datapath, regfile, nustar_path, edit_regfile=edit_regfile, 
                                      nofit=nofit, clobber=True, pile_up_corr=pile_up_corr, nuradius=nuradius)
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
                e_lo1, e_hi1, eff_area = nuutil.read_arf(arf_files[0])
                e_lo2, e_hi2, rmf_mat = nuutil.read_rmf(rmf_files[0])
                engs,cnts,lvtm,ontim=nuutil.read_pha(pha_files[0])
                if pile_up_corr:
                    arf_files_unphys, rmf_files_unphys, pha_files_unphys = find_nuproducts(nustar_path, timestring, fpm,
                                                                               special_pha=special_pha, grade='21_24')
                    e_lo1, e_hi1, eff_area_u = nuutil.read_arf(arf_files_unphys[0])
                    e_lo2, e_hi2, rmf_mat_u = nuutil.read_rmf(rmf_files_unphys[0])
                    engs,cnts_u,lvtm,ontim=nuutil.read_pha(pha_files_unphys[0]) 
                    
                    
        else:
            print('You need to set make_nustar=True to make spectral data products.')
            return
    
    print('')
    
    #Get area of NuSTAR region
    if nofit:
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
            
            cnts_corr = cnts-(5./4)*cnts_u
        else:
            plt.plot(engs, (cnts-0.25*cnts_u), label='corr')
            plt.plot(engs, cnts, label='og grade 0')
            plt.plot(engs, cnts_u, label='unphys')
            plt.yscale('log')
            plt.xlim([0,20])
            plt.legend()
            plt.savefig(nustar_path+timestring+'/'+timestring+fpm+'pile_up.png')
        
            cnts_corr = cnts-0.25*cnts_u
    
    numax = enz[-1]
    
    print('Max NuSTAR Energy: ', numax)
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
    
    if actual_total_counts:
        return rate, erate, nutrs, tresp, logt, fpm, atc
    else:
        return rate, erate, nutrs, tresp, logt, fpm


def edit_nuscreen(path_to_dodem, nustar_path, timestring, fpm, datapath):
    """
    Requires shell script run_nuscreen.sh in input path (nustar_path). 
    Edits based on inputs + overwrites new version.
    """
    
    with open(path_to_dodem+'run_nuscreen.sh', "r+") as f:
        screen = f.read()
        #print(screen)

        screenlist = screen.split('\n')
        #print(screenlist)
        screenlist[4] = 'interval='+timestring
        screenlist[5] = 'fpm='+fpm
        screenlist[6] = 'INDIR='+datapath
        screenlist[7] = 'work_dir='+nustar_path

        newscreen = '\n'.join(screenlist)
        f.seek(0)
        f.write(newscreen)
        f.truncate()
            
            
def edit_nuproducts(path_to_dodem, nustar_path, timestring, fpm, regfile, datapath, unphys_products=0):
    """
    Requires shell script run_nuproducts.sh in input path (nustar_path). 
    Edits based on inputs + overwrites new version.
    """
    
    with open(path_to_dodem+'run_nuproducts.sh', "r+") as f:
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
            print(outfile, 'exists! Not re-making it.')
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


def find_intervals(macro_interval, nuenergies, datapath, filename, gtifile='', 
                   regfile='', nofit=False, firstrun=False):
    """
    For a given uninterupted NuSTAR data interval (say, an orbit), find a list of time intervals for 
    DEMs with duration ranging from 30s to however long is needed to contain 10 actual NuSTAR counts 
    in a supplied energy range.
    
    Keywords
    ---------
    
    macro_interval - full data interval (i.e. orbit) (tuple of astropy Time objects)
                FORMAT LIKE, 
                time=(astropy.time.Time('2018-05-29T19:08:00', scale='utc'), 
                        astropy.time.Time('2018-05-29T19:14:00', scale='utc')
    
    nuenergies - energy range to check for sufficient statistics
                  FORMAT LIKE,
                  nuenergies = [6.,10.]
                  
    filename - component of name of picklefile where the optimized time intervals will be saved.
                  
    gtifile - template file, will be changed based on time range.
    
    datapath - path to NuSTAR data directory (top level, i.e. obsid-named directory)
    
    regfile - region file. Will be edited (optimal region found) if edit_regfile==True.
    
    nofit - Set nofit=True to use center of mass regions (not caring about chip gap), False
    firstrun - Set firstrun=True to clobber existing spectral data products.
    
    """
    
    if gtifile=='':
        print('We require a starter GTI file to edit (can be found in datapath+/event_cl).')
        return
    if regfile=='':
        print('We require a starter circular region file to edit (can be made manually in ds9).')
        return
        
              
    nustar_path='./'
    baseinterval = 30*u.s
    time_intervals=[]
    time1 = macro_interval[0]
    time2 = macro_interval[0]+baseinterval
    
    
    while time2 <= macro_interval[1]:
        time_int = [time1, time2]
        print(time_int)
        res = combine_fpm(time_int, nuenergies, nustar_path, make_nustar=True,
                                              gtifile=gtifile, datapath=datapath, regfile=regfile, 
                                                edit_regfile=True, actual_total_counts=True,
                                                nofit=nofit, clobber=firstrun)
        if res[1]:
            print('Enough Counts!')
            time_intervals.append(time_int)
            print(time_int)
            time1=copy.deepcopy(time2)
            time2+=baseinterval
        else:
            counts=res[0]
            print('only', counts)
            if counts > 5:
                print('Extending 15s...')
                time2+=baseinterval/2
            else:
                print('Extending 30s...')
                time2+=baseinterval


    #NOT DONE HANDLING FINAL INTERVAL
    #Time2 has just increased, and crossed the end of the macro interval. This means either that we've 
    #just sucessfully completed a full interval and are starting a new one, or that we've extended time2
    #trying to complete an interval in=progress.

    #New last interval: from time1 to the end
    last_interval = [time1, macro_interval[1]]
    lastres = combine_fpm(last_interval, nuenergies, nustar_path, make_nustar=True,
                                                gtifile=gtifile, datapath=datapath, regfile=regfile,
                                                edit_regfile=True, actual_total_counts=True,
                                                nofit=nofit, clobber=firstrun)
    #If the new last interval is complete, add it.
    if lastres[1]:
        time_intervals.append(last_interval)
    #If not, edit the current last interval to extend its end to the end of the macro interval.
    else:
        time_intervals[-1][1] = macro_interval[1]


    for t in time_intervals:
        print(t[0].strftime('%H-%M-%S'), t[1].strftime('%H-%M-%S'))
        
    data = {'time_intervals': time_intervals}
    
    if nofit:
        filename='noregionfit_'+filename

    with open(filename, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)       

    return time_intervals




def get_saved_intervals(nofit=False, basedir='./', custom_file=[]):
    """
    Reads in a file of the type made by find_intervals() above, containing DEM time intervals.
    
    Keywords
    ---------
    
    nofit - If True, checks for filename of the type find_intervals() gives when NOT using region fitting.
            If False, checks for filename of the type find_intervals() gives when using region fitting.
            
    basedir - path to where the file is located. 
    
    custom_file - Set to a specific file to use (still within the basedir).
    
    """
    
    
    if nofit==True:
        time_interval_files = glob.glob(basedir+'/noregionfit_time_intervals*.pickle')
    else:
        time_interval_files = glob.glob(basedir+'/time_intervals*.pickle')
        
    if bool(custom_file):
        time_interval_files = glob.glob(basedir+custom_file)
        
    print(len(time_interval_files))
    
    #print(time_interval_files)
    num=[]
    for i in range(0,len(time_interval_files)):
        filename = time_interval_files[i].split('/')[-1]
        nunum = [int(s) for s in [*filename] if s.isdigit()]
        if len(nunum) > 1:
            print('Not set up to sort file names with more than one number!')
            print('Culprit: ', filename)
            return
        num.extend(nunum)
    sort_index = np.argsort(num)
    time_interval_files = [time_interval_files[s] for s in sort_index] 
                          
    #print(time_interval_files)
    time_intervals=[]
    for tf in time_interval_files:
        with open(tf, 'rb') as f:
            data = pickle.load(f)
            if type(data) == tuple:
                time_intervals.extend(data[0]['time_intervals'])
            if type(data) == dict:
                time_intervals.extend(data['time_intervals'])
            
    return time_intervals



