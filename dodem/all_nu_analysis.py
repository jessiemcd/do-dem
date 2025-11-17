import numpy as np
from astropy import units as u
import pickle
import pathlib

path_to_dodem = '/Users/jmdunca2/do-dem/'
from sys import path as sys_path
sys_path.append(path_to_dodem+'/dodem/')

import time_interval_selection as tis
import visualize_dem_results as viz
import gauss2D as g2d
import images_and_coalignment as iac
import nustar_dem_prep as nu


def get_nustar_orbit_times(datapath):

    import nustar_utilities as nuutil

    evt_data, hdr = nu.return_submap(datapath=datapath, fpm='A', return_evt_hdr=True)
    time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
    
    return time0, time1


def goes_string(flux):

    #goes_conversion: [1e-8, 1e-7, 1e-6, 1e-5, 1e-4] = ['A', 'B', 'C', 'M', 'X']
    goes_levels = np.array([-12,-11,-10, -9, -8, -7, -6, -5, -4])
    goes_labels = np.array(['A0.0001', 'A0.001', 'A0.01','A0.1', 'A', 'B', 'C', 'M', 'X'])

    order_ = omag(flux.value)
    ind = np.where(goes_levels == order_)[0]
    if order_ < -8:
        round_ = -1*(order_+7)
        numberstr = str(round(flux.value/(10**(-8)), round_)) 
        gstring = 'A'+numberstr
    else:
        gstring = goes_labels[ind][0]+str(round(flux.value/(10**(order_)), 1)) 

    return gstring

def goes_minmax(vals):
    """
    Takes in a list of GOES XRSB values and returns the minimum and maximum values, 
    both as numbers and nicely formatted strings.
    
    """
    gmin, gmax = np.nanmin(vals), np.nanmax(vals)

    if gmin < 0:
        print('Negative?', gmin)
    if gmax < 0:
        print('Negative?', gmax)

    
    goesrange_vals = [gmin.value, gmax.value]*gmin.unit
    goesrange_strings = [goes_string(gmin),
                         goes_string(gmax)]
    
    print('GOES Class range: ', goesrange_strings[0], \
          ' to ', goesrange_strings[1])

    return goesrange_vals, goesrange_strings


def try_get_goes(datapath, satellite=13):
    """
    Get GOES data + return minimum and maximum GOES class during a NuSTAR orbit
    """

    #GOES16 data seems to start in february 2017 
    
    #Path to top-level do-dem directory - edit for your system.
    path_to_dodem = '/Users/jmdunca2/do-dem/'
    from sys import path as sys_path
    sys_path.append(path_to_dodem+'/dodem/')

    import nustar_utilities as nuutil
    import lightcurves as lc
    import astropy.time

    import datetime
    
    fpm='A'
    #Get first and last times associated with an event in the cleaned event file.
    evt_data, hdr = nu.return_submap(datapath=datapath, fpm=fpm, return_evt_hdr=True)
    time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]

    time0_ = time0.tt.datetime
    time1_ = time1.tt.datetime

    #print(time0, time1)

    timerange = [time0_,time1_]

    #from datetime import timezone
    #timerange = [t.replace(tzinfo=timezone.utc) for t in timerange]

    lc.get_goes(timerange, satellite=satellite, peek=False)
    instrument='GOES'
    data = lc.load_lightcurves(instrument)
    
    ylabel = data['GOES flux label']
    goestimes = data['GOES Times']
    xrsbcounts = data['XRSB counts']
    xrsblabel = data['XRSB label']
    gts = np.array([t.datetime for t in goestimes])

    int_inds = np.where(np.logical_and(gts > time0_, gts < time1_))[0]
    int_counts = xrsbcounts[int_inds]

    return goes_minmax(int_counts)


def all_obs_goes(datapaths, satellite=13):
    """
    Compile GOES info for all orbits in a given AR observation.
    """

    all_goes_vals=[]
    goes_per_orbit=[]
    goes_per_orbit_strings=[]
    for id in datapaths:
        goesrange_vals, goesrange_strings = try_get_goes(id, satellite=satellite)
        all_goes_vals.extend(goesrange_vals)
        goes_per_orbit.append(goesrange_vals)
        goes_per_orbit_strings.append(goesrange_strings)


    print('')
    print('all-day:')
    all_goes_vals_ = [vv.value for vv in all_goes_vals]*all_goes_vals[0].unit
    allminmax, allstrings = goes_minmax(all_goes_vals_)

    return allminmax, allstrings, goes_per_orbit, goes_per_orbit_strings

    
    

def omag(x):
    return int(np.floor(np.log10(np.abs(x))))


def get_durations(datapaths, fpm='B', filepaths=False, hkpaths=[]):

    """
    If filepaths == False, expects datapaths to be links to NuSTAR data directories.

    If filepaths == True, expects datapaths to lead to SPECIFIC SUNPOS EVT FILES, and hkpaths 
        must be added to host the NuSTAR data directories. They must be corresponding
        lists of the same length.
    """


    #Path to top-level do-dem directory - edit for your system.
    path_to_dodem = '/Users/jmdunca2/do-dem/'
    from sys import path as sys_path
    sys_path.append(path_to_dodem+'/dodem/')
    
    import nustar_utilities as nuutil

    from astropy import units as u
    import glob
    from astropy.io import fits

    import astropy.time
    


    durations = []
    lvttotals = []
    total = 0.*u.min
    for d in range(0, len(datapaths)):
        id = datapaths[d]
        #print(id)
        #Get first and last times associated with an event in the cleaned event file.
        if filepaths:
            evt_data, hdr = nu.return_submap(specific_evt=id, fpm=fpm, return_evt_hdr=True, 
                                             already_sunpos=True)
            hk = glob.glob(hkpaths[d]+'/hk/*'+fpm+'_fpm.hk')
        else:
            evt_data, hdr = nu.return_submap(datapath=id, fpm=fpm, return_evt_hdr=True)
            hk = glob.glob(id+'/hk/*'+fpm+'_fpm.hk')
            
        time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
        #Duration: difference between them.
        durations.append((time1-time0).to(u.min))

        #Load in housekeeping data
        hdulist = fits.open(hk[0])
        dat = hdulist[1].data
        hdr = hdulist[1].header
        hdulist.close()

        #Get the livetimes from between the time0,1 found above.
        livetimes = dat['livetime']
        mjd_ref_time=astropy.time.Time(hdr['mjdrefi'],format='mjd')
        livetime_times = astropy.time.Time(mjd_ref_time+dat['time']*u.s,format='mjd')
        durbins = np.where(np.logical_and(livetime_times > time0, livetime_times < time1))[0]
        livetimes_ = np.array(livetimes)[durbins]
        
        #print(id[-12:-1], (np.sum(livetimes_)*u.s)/((time1-time0).to(u.s)))

        lvttotals.append((np.sum(livetimes_)*u.s).to(u.min))

    return durations, sum(durations), sum(lvttotals)


def print_times(files):

    #Path to top-level do-dem directory - edit for your system.
    path_to_dodem = '/Users/jmdunca2/do-dem/'
    from sys import path as sys_path
    sys_path.append(path_to_dodem+'/dodem/')
    
    import nustar_utilities as nuutil
    
    for d in files:
        evt_data, hdr = nu.return_submap(datapath=d, fpm='A', return_evt_hdr=True)
        time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
        print(time0, time1)



def load_nufiles(f):
    """
    """
    hdulist = fits.open(f)
    dat = hdulist[1].data
    hdr = hdulist[1].header
    hdulist.close()
        
    return dat, hdr


def get_exposures(target_dict, dogoes=False):

    keys = target_dict.keys()
    
    all_effective_exposure=0*u.min
    all_duration=0*u.min
    all_all_goes=[]
    for key in keys:
        ARdict = target_dict[key]
        durations, total, lvttotal = get_durations(ARdict['datapaths'])
        print(key)
        print('Duration: ', np.round(total, 2), '| Effective Exposure: ', np.round(lvttotal,2), '| Livetime %: ', np.round((lvttotal/total),2))
    
        target_dict[key]['orbit durations'] = durations
        target_dict[key]['total duration'] = total
        target_dict[key]['total livetime'] = lvttotal
    
        all_effective_exposure += lvttotal
        all_duration += total
    
        if dogoes:
            allminmax, allstrings, \
                goes_per_orbit, goes_per_orbit_strings = all_obs_goes(ARdict['datapaths'], satellite=ARdict['goes_satellite'])
            target_dict[key]['AR GOES min, max vals'] = allminmax
            target_dict[key]['AR GOES min, max strings'] = allstrings
            target_dict[key]['orbit GOES min, max vals'] = goes_per_orbit
            target_dict[key]['orbit GOES min, max strings'] = goes_per_orbit_strings
            #print('')
            all_all_goes.append(allminmax)
    
    
    #print('')
    print('Total Dataset Effective Exposure: ', all_effective_exposure.to(u.h))
    print('Total Dataset Observation Duration: ', all_duration.to(u.h))
    
    if dogoes:
        print('GOES - all-dataset: ')
        all_all_goes_vals_ = [vv.value for vv in all_all_goes]*all_all_goes[0].unit
        allminmax, allstrings = goes_minmax(all_all_goes_vals_)




def single_gauss_prep(key, file, plot=True, guess=[], make_scripts=True,
                     plotregion=[], plotgaussregions=False):


    with open(file, 'rb') as f:
        data = pickle.load(f)
    
    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    working_dir = ARDict['working_dir']
    
    #Make a new working directory for prepped data/etc if it doesn't yet exist
    save_path = pathlib.Path(working_dir)
    if not save_path.exists():
        save_path.mkdir()

    gauss_stats=[]
    for i in range(0, len(id_dirs)):
        #guess, fast_min_factor 
        res = g2d.per_orbit_onegauss_params(id_dirs[i], guess=guess, plot=plot,
                                           plotregion=plotregion, plotgaussregions=plotgaussregions)
        gauss_stats.append(res)


    ARDict['gauss_stats'] = gauss_stats

    data[key] = ARDict
    
    with open(file, 'wb') as f:
             # Pickle the 'data' dictionary using the highest protocol available.
             pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 

    if make_scripts:
        #where: where to find templates + place scripts.
        tis.make_tis_scripts(obsids, key, where='./scripts/')  



def double_gauss_prep(key, file, plot=True, guess=[], guess2=[], sep_axis='SN', make_scripts=True,
                      plotregion=[], write_input_regions=False,
                      plotgaussregions=False, write_regions=False):


    with open(file, 'rb') as f:
        data = pickle.load(f)
    
    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    working_dir = ARDict['working_dir']
    
    #Make a new working directory for prepped data/etc if it doesn't yet exist
    save_path = pathlib.Path(working_dir)
    if not save_path.exists():
        save_path.mkdir()

    gauss_stats=[]
    for i in range(0, len(id_dirs)):
        #guess, fast_min_factor 
        res = g2d.per_orbit_twogauss_params(id_dirs[i], sep_axis=sep_axis, guess=guess, guess2=guess2, plot=plot,
                                           plotregion=plotregion, plotgaussregions=plotgaussregions,
                                            write_input_regions=write_input_regions,
                                            write_regions=write_regions, region_dir=working_dir)
                        
        gauss_stats.append(res)
        #print('')

    print(gauss_stats)


    ARDict['gauss_stats'] = gauss_stats

    data[key] = ARDict
    
    with open(file, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 

    if make_scripts:
        ##where: where to find templates + place scripts.
        tis.make_tis_scripts(obsids, key, where='./scripts/', tworegion=True)   
    



def manual_prep(key, file, plot=True, make_scripts=True,
                      inputregion=[], write_input_regions=True,
                      plotgaussregions=False):

    """
    key - key for all target dictionary (where to get information about the nustar data, region, etc)
    file - file containing all target directory
    plot - set True for plots to be made in general

    make_scripts - set True with you've finalized your regions and are ready to write corresponding scripts for TIS
    inputregion - list of region dictionaries, of the form:
                    
                    plotregion = [{'centerx': 950, 'centery': -325, 'radius': 150},
                               {'centerx': 900, 'centery': -50, 'radius': 150}]

                (You can have as many regions as you want. Values in arcseconds from solar center.)
                
    write_input_regions - set True to write a .reg file for every region in plotregion. This is needed to run TIS 
                            (scripts written in make_scripts will call functions that will look for these regions).

    plotgaussregions - set True to plot 150" circles centered at gaussian fit result locations, if useful for visualization.

    
    """
    guess=[]
    guess2=[]


    with open(file, 'rb') as f:
        data = pickle.load(f)
    
    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    working_dir = ARDict['working_dir']

    
    #Make a new working directory for prepped data/etc if it doesn't yet exist
    save_path = pathlib.Path(working_dir)
    if not save_path.exists():
        save_path.mkdir()

    region_stats=[]
    for i in range(0, len(id_dirs)):
        #guess, fast_min_factor 
        res = g2d.per_orbit_manual_params(id_dirs[i], guess=guess, guess2=guess2, plot=plot,
                                           plotregion=inputregion, plotgaussregions=plotgaussregions,
                                            write_input_regions=write_input_regions,
                                            region_dir=working_dir)
                        
        region_stats.append(res)
        #print('')


    ARDict['region_stats'] = region_stats

    data[key] = ARDict
    
    with open(file, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 

    if make_scripts:
        ##where: where to find templates + place scripts.
        tis.make_tis_scripts(obsids, key, where='./scripts/', manualregion=True)  
    



def do_key_dem(key, file, 
               missing_last=False, missing_orbit=4, 
               plot_xrt=True, plot_aia=False,
               use_prepped_aia=True,
               high_temp_analysis=False, rscl=True, 
               do_no_xrt_version=False,
               do_onlyaia_version=False,
               do_aiaxrt_version=False,
              aia_region_dict={}, input_time_intervals=[],
              pick_region=False, regionind=0,
              save_inputs_file=False):

    """
    Keywords:
    -----------
    
    key - obs. key (indexes dictionary contained in file)
    file - saved by AR_inventory, contains dictionary of information about each included observation
    
    high_temp_analysis - set to do multiple DEMs, varying the temperature bounds at the high side
    rscl - set True to do the rscl process
    
    missing_last - Set missing_last=True to trim time interval list to exclude the last interval in an orbit (missing_orbit)
    (useful due to some NCCS AIA data intervals slightly shorter than NuSTAR intervals). 
    missing_orbit - the orbit for which to apply missing_last
    
    plot_aia, plot_xrt - passed to dodem, make plots of data for reference/troubleshooting.
    do_no_xrt_version, do_onlyaia_version, do_aiaxrt_version - set each True to do this kind of DEM. Only one should be set, or else they
                                                                override eachother in order. 
                                                                
    save_inputs_file - Set to save DEM inputs, but NOT do DEM. Only works with high_temp_analysis=False (direct dodem call)
    
    AIA OPTIONS
    use_prepped_aia - set True to use the AIA inputs saved in the file listed in the information dictionary for this key under 
                        the key 'prepped_aia'.
                       This is typically the method used with NCCS-output AIA files.
     * if ^ False *
    aia_region_dict - set to dictionary of AIA regions to download + prep AIA data while doing dem (new cutout implementation)
    pick_region - set True to choose only one region to do DEMs of, within the above list.
    regionind - the index of the region chosen (only used if pick_region = True)
    
    input_time_intervals - set to list of chosen time intervals to override the list of all available time intervals found via TIS.
    
    
    
    Other Notes:
    
    
    """


    #Path to top-level do-dem directory - edit for your system.
    path_to_dodem = '/Users/jessieduncan/do-dem/'
    
    #Empty, so function call later doesn't fail even if these aren't updated. 
    input_aia_region=""
    input_aia_region_dict={}
    fetch_cutout=False

    with open(file, 'rb') as f:
        data = pickle.load(f)

    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    working_dir = ARDict['working_dir']
    if use_prepped_aia:
        prepped_aia_dir = ARDict['prepped_aia']
    method=ARDict['method']
    if method=='double':
        gauss_stats = ARDict['gauss_stats']
        sep_axis = gauss_stats[0][0]
    else:
        sep_axis = ''
        


    if method in ['input', 'double']:
        directories = get_region_directories(key, targets_file=file)
        print('1',directories)
        all_all_time_intervals, fixit = tis.region_time_intervals(directories, id_dirs, shush=True)        
        if pick_region:
            directories = [directories[regionind]]
            print('2', directories)

        

    if method=='fit':
        onegauss=True
        regfile=path_to_dodem+'starter_region.reg'
        all_time_intervals, all_time_intervals_list = tis.find_all_intervals(working_dir, shush=True, 
                                                                        missing_last=missing_last, missing_orbit=missing_orbit)

    
    if input_time_intervals:
        if method in ['input', 'double']:
             all_all_time_intervals = input_time_intervals
        
        if method=='fit':  
            all_time_intervals = input_time_intervals
        
    
    #What instruments are you using?
    #---------------------------------
    aia=True
    #---------------------------------
    eis=False
    #This is where I'm putting my XRT level-1 data and grade maps:
    xrt_path=working_dir+'/XRT_for_DEM/'
    ogxrt=True
    
    from astropy import units as u
    exposure_dict={'Be_thin': [1*u.s, 100*u.s],
                    'Be_med': [],
                  'Al_poly': [0.1*u.s, 100*u.s]}
    #---------------------------------
    plot=False
    #---------------------------------
    nustar=True
    force_nustar=True
    #force_nustar=False
    #nustar=False
    
    #If nustar is being used, here are the chosen energy ranges:
    nuenergies=[[2.5,3.5], [3.5,6.], [6.,10.]]
    nuradius=150
    #---------------------------------
    
    #---------------------------------
    #---------------------------------
    #What temperature range would you like to use? (units: log(T))
    minT=5.6
    maxT=7.2
    
    #---------------------------------
    #---------------------------------

    name=key
    
    if len(np.where([do_no_xrt_version, do_onlyaia_version, do_aiaxrt_version])[0]) > 1:
        print('Please choose only one of do_no_xrt_version, do_onlyaia_version, do_aiaxrt_version.')
        return

    if do_no_xrt_version:
        name=key+'_no_xrt'
        ogxrt=False
        
    if do_onlyaia_version:
        name=key+'_onlyaia'
        ogxrt=False
        force_nustar=False
        nustar=False

    if do_aiaxrt_version:
        name=key+'_aiaxrt'
        force_nustar=False
        nustar=False
        
    import dodem
    import glob
    
    for o in range(0, len(obsids)):
        xrt=ogxrt
        datapath=id_dirs[o]
        xrt_path=path_to_dodem+'other_idl/'+obsids[o]+'_coobs/XRT_for_DEM/'
        print(xrt_path)
        if not pathlib.Path(xrt_path).is_dir():
            xrt=False
        gtifile=datapath+'event_cl/nu'+obsids[o]+'A06_gti.fits'
        if use_prepped_aia:
            orbit_aia_dir = prepped_aia_dir+'/orbit_'+obsids[o]+'/'
        obsid=obsids[o]

        if method=='fit':
            guess = ARDict['gauss_stats'][o][0]
            time_intervals = all_time_intervals[o] 
            #If downloading new AIA, put it in DEM folder.
            sunpy_dir=working_dir
            for time in time_intervals:

                if use_prepped_aia:
                    #print(orbit_aia_dir)
                    sunpy_dir=working_dir
                    res = iac.read_interval_dicts(time, where=orbit_aia_dir, bltr=True)
                    if res is None:
                        print('Found no AIA')
                        continue
                    else:
                        data, bl, tr, region_input = res
                        #print(region_input)
                        if 'region0' in data.keys():
                            #print(data.keys())
                            data = data['region0']
                            region_input = region_input[0]
                        #print(data['aia_dn_s_px'])
                        input_aia_region=[]
                        input_aia_region_dict=[]

                elif aia_region_dict:
                    fetch_cutout=True
                    data=[]
                    input_aia_region="circle"
                    input_aia_region_dict={'center': [aia_region_dict[o]['centerx'], aia_region_dict[o]['centery']],
                                           'radius': aia_region_dict[o]['radius']*u.arcsec
                                            }
                    region_input=input_aia_region_dict
                    bl, tr = [aia_region_dict[o]['centerx']-(300*u.arcsec), aia_region_dict[o]['centery']-(300*u.arcsec)], \
                            [aia_region_dict[o]['centerx']+(300*u.arcsec), aia_region_dict[o]['centery']+(300*u.arcsec)]
                else:
                    print('Did not specify a working AIA method.')
                    bl, tr = [], []
                    data = []
                    region_input=[]
                    aia=False
                    
                
                if high_temp_analysis:
                    dodem.high_temp_analysis(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, name2=name,
                                                   highT=7.2, #(minT, maxT are varied)
                                                   working_directory=working_dir, #(plotresp set false in high_temp_analysis)
                                                   default_err=0.2, path_to_dodem=path_to_dodem,
                                                   demmethod='DEMREG', use_prior_prep=False,
                            
                                                   #demreg/xrt_iterative related
                                                   rgt_fact=1.2, max_iter=30, rscl=rscl,
                                                   reg_tweak=1, #mc_in, mc_rounds hard coded in high_temp_analysis (same as below)
                            
                                                   #nustar=related
                                                   nuenergies = nuenergies, #combine_fpm hard coded in high_temp_analysis (same as below)
                                                   datapath=datapath, gtifile=gtifile, 
                                                   nuradius=nuradius, guess=guess, onegauss=onegauss,
                                                   adjacent_grades=True, pile_up_corr=True,
                                                   force_nustar=force_nustar,
                                             
                                                   #aia related
                                                   load_prepped_aia=data, sunpy_dir=sunpy_dir,
                                                   input_aia_region=input_aia_region, input_aia_region_dict=input_aia_region_dict,
                                                   real_aia_err=True, fetch_cutout=fetch_cutout,
                            
                                                   #xrt related
                                                   xrtmethod='Average', real_xrt_err=True, xrt_path=xrt_path,
                                                   xrt_exposure_dict = exposure_dict, plot=plot_xrt, #(this plots xrt)
                                                   input_xrt_region="circle", input_xrt_region_dict=region_input)
                    


                else:
                    #print('xrt is: ', xrt)
                    dodem.dodem(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, name=name,
                                            minT=minT, maxT=maxT,
                                            working_directory=working_dir,
                                            default_err=0.2, path_to_dodem=path_to_dodem,
                                            save_inputs_file=save_inputs_file,
                    
                                            #demreg related
                                            rgt_fact=1.2, max_iter=30, rscl=rscl,
                                            reg_tweak=1, mc_in=True, mc_rounds=100, 
                                            
                                            #nustar related 
                                            combine_fpm=True, nuenergies=nuenergies, 
                                            datapath=datapath, gtifile=gtifile,
                                            nuradius=nuradius, guess=guess, onegauss=onegauss,
                                            adjacent_grades=True, pile_up_corr=True,
                                            force_nustar=force_nustar,
                    
                                            #aia related
                                            load_prepped_aia=data, sunpy_dir=sunpy_dir,
                                            input_aia_region=input_aia_region, input_aia_region_dict=input_aia_region_dict,
                                            real_aia_err=True, fetch_cutout=fetch_cutout, plot_aia=plot_aia,
        
                                            #xrt related
                                            xrtmethod='Average', real_xrt_err=True, xrt_path=xrt_path,
                                            xrt_exposure_dict=exposure_dict, plot_xrt=plot_xrt,
                                            input_xrt_region="circle", input_xrt_region_dict=region_input)
        #print('')

        if method in ['input', 'double']:
            print('xrt is: ', xrt)
            fpm='A'
            if method=='input':
                regfiles = glob.glob(working_dir+'gauss_cen_'+obsid+'_'+fpm+'_user_input*.reg')
            if method=='double':
                regfiles = glob.glob(working_dir+'gauss_cen_'+obsid+'_'+fpm+'_*.reg')
                
            regfiles.sort()
            if pick_region:
                regfiles = [regfiles[regionind]]
            
            for i in range(0, len(directories)):
                #Time intervals for this region, orbit
                time_intervals = all_all_time_intervals[i][o]

                regfile = regfiles[i]
                print(directories[i], regfiles[i])
                #If downloading new AIA, put it in DEM folder.
                sunpy_dir = directories[i]


                for time in time_intervals:
                    if use_prepped_aia:
                        res = iac.read_interval_dicts(time, where=orbit_aia_dir, bltr=True)
                        if res is None:
                            print('Found no AIA')
                            continue
                        datas, bl, tr, xrt_region_inputs = res
                        data = datas['region'+str(i)]
                        region_input = xrt_region_inputs[i]
                        
                    elif aia_region_dict:
                        fetch_cutout=True
                        data=[]
                        aia_region_dict_=aia_region_dict[o][i]
                        input_aia_region="circle"
                        input_aia_region_dict={'center': [aia_region_dict_['centerx'], aia_region_dict_['centery']],
                                               'radius': aia_region_dict_['radius']*u.arcsec
                                                }
                        region_input=input_aia_region_dict
                        apdw = 300
                        bl, tr = [aia_region_dict_['centerx']-(apdw*u.arcsec), aia_region_dict_['centery']-(apdw*u.arcsec)], \
                                [aia_region_dict_['centerx']+(apdw*u.arcsec), aia_region_dict_['centery']+(apdw*u.arcsec)]
                    else:
                        print('Did not specify a working AIA method.')
                        bl, tr = [], []
                        data = []
                        region_input=[]
                        aia=False



                    if high_temp_analysis:
                        dodem.high_temp_analysis(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, name2=name,
                                                highT=7.2, #(minT, maxT are varied)
                                                working_directory=directories[i], #(plotresp set false in high_temp_analysis)
                                                default_err=0.2, path_to_dodem=path_to_dodem,
                                                demmethod='DEMREG', use_prior_prep=False,

                                 
                                                #demreg/xrt_iterative related
                                                rgt_fact=1.2, max_iter=30, rscl=rscl,
                                                reg_tweak=1, #mc_in, mc_rounds hard coded in high_temp_analysis (same as below)

                            
                                                #nustar=related
                                                nuenergies = nuenergies, #combine_fpm hard coded in high_temp_analysis (same as below)
                                                datapath=datapath, gtifile=gtifile, 
                                                nuradius=nuradius, edit_regfile=False,
                                                regfile=regfile,
                                                adjacent_grades=True, pile_up_corr=True,
                                                force_nustar=force_nustar,

                                                #aia related
                                                load_prepped_aia=data, sunpy_dir=sunpy_dir,
                                                input_aia_region=input_aia_region, input_aia_region_dict=input_aia_region_dict,
                                                real_aia_err=True, fetch_cutout=fetch_cutout,

                                                #xrt related
                                                xrtmethod='Average', real_xrt_err=True, xrt_path=xrt_path,
                                                xrt_exposure_dict=exposure_dict, plot=plot_xrt, #(this plots xrt)
                                                input_xrt_region="circle", input_xrt_region_dict=region_input)

                                                 
                    else:
                        print('xrt is: ', xrt)
                        dodem.dodem(time, bl, tr, xrt=xrt, aia=aia, nustar=nustar, name=name,
                                    minT=minT, maxT=maxT,
                                    working_directory=directories[i],
                                    default_err=0.2, path_to_dodem=path_to_dodem,
                                    save_inputs_file=save_inputs_file,
            
                                    #demreg related
                                    rgt_fact=1.2, max_iter=30, rscl=rscl,
                                    reg_tweak=1, mc_in=True, mc_rounds=100, 
                                    
                                    #nustar related 
                                    combine_fpm=True, nuenergies=nuenergies, 
                                    datapath=datapath, gtifile=gtifile,
                                    nuradius=nuradius, edit_regfile=False,
                                    regfile=regfile,
                                    adjacent_grades=True, pile_up_corr=True,
                                    force_nustar=force_nustar,
            
                                    #aia related
                                    load_prepped_aia=data, sunpy_dir=sunpy_dir,
                                    input_aia_region=input_aia_region, input_aia_region_dict=input_aia_region_dict,
                                    real_aia_err=True, fetch_cutout=fetch_cutout, plot_aia=plot_aia,
        
                                    #xrt related
                                   xrtmethod='Average', real_xrt_err=True, xrt_path=xrt_path,
                                    xrt_exposure_dict=exposure_dict, plot_xrt=plot_xrt,
                                    input_xrt_region="circle", input_xrt_region_dict=region_input)



def get_key_resultfiles(key, file, fromhome=False,
                        withparams=False,
                       namesearchstring='',
                       shush=False):


    
    import glob
    #print('doing ', key)

    if fromhome:
        
        with open(file, 'rb') as f:
            data = pickle.load(f)

        ARDict = data[key]
        
        id_dirs = ARDict['datapaths']
        #obsids = ARDict['obsids']
        working_dir = ARDict['working_dir']
        # prepped_aia_dir = ARDict['prepped_aia']
        method = ARDict['method']
        
        if method in ['input', 'double']:
            directories = get_region_directories(key, targets_file=file)
            #all_all_time_intervals is a list... 
            #       with entries (lists) for each directory/region
            #             those lists have entries (lists) for each orbit
            #                   those lists have entries for each time interval.
            all_all_time_intervals, fixit = tis.region_time_intervals(directories, id_dirs, shush=True)



            res_files=[]
            for d in range(0, len(directories)):
                d_files=[]
                dir_ = directories[d]
                all_orbits_time_intervals = all_all_time_intervals[d]
                for ot in range(0, len(all_orbits_time_intervals)):
                    orbittimes = all_orbits_time_intervals[ot] 
                    if not orbittimes:
                        continue
                        
                    for tt in orbittimes:
                        timestring = viz.make_timestring(tt)
                        if withparams:
                            files = glob.glob(dir_+'/'+timestring+'/'+'*5.6_7.2*'+namesearchstring+'*withparams.pickle')
                        else:
                            #print(dir_+'/'+timestring+'/'+'*5.6_7.2*'+namesearchstring+'*.pickle')
                            fs = glob.glob(dir_+'/'+timestring+'/'+'*5.6_7.2*'+namesearchstring+'*.pickle')
                            files = [f for f in fs if 'withparams' not in f]
                            
                        try:
                            d_files.append(files[0])
                        except IndexError:
                            if not shush:
                                print(key, d, timestring, files)
                        
                if len(directories) > 1:
                    res_files.append(d_files)
                else:
                    return d_files, False

            return res_files, True

        
        if method=='fit':
            res_files=[]
            all_time_intervals, all_time_intervals_list = tis.find_all_intervals(working_dir, shush=True)
            for ot in range(0, len(all_time_intervals)):
                orbittimes = all_time_intervals[ot]  
                for tt in orbittimes:
                    timestring = viz.make_timestring(tt)
                    if withparams:
                        files = glob.glob(working_dir+'/'+timestring+'/'+'*5.6_7.2*'+namesearchstring+'*withparams.pickle')
                    else:
                        fs = glob.glob(working_dir+'/'+timestring+'/'+'*5.6_7.2*'+namesearchstring+'*.pickle')
                        files = [f for f in fs if 'withparams' not in f]
                        #print(files)
                        #print('')
                        
                    try:
                        res_files.append(files[0])
                    except IndexError:
                        pass
                        #print(key, timestring, files)
                        
            return res_files, False
            
    else:
        if withparams:
            res_files = glob.glob('./compact_results/*'+key+'*'+namesearchstring+'*withparams.pickle')
        else:
            rfs = glob.glob('./compact_results/*'+key+'*'+namesearchstring+'*.pickle')
            res_files = [f for f in rfs if 'withparams' not in f]
            #print(res_files)
            
        res_files.sort()
            
        if 'region' in res_files[0]:
            rez = []
            zeros = [f for f in res_files if 'region_0' in f]
            if zeros:
                rez.append(zeros)
            ones = [f for f in res_files if 'region_1' in f]
            if ones:
                rez.append(ones)
    
    
            if zeros and ones:
                return rez, True
            else:
                #print(rez[0])
                return rez[0], False
    
    
        else:
            return res_files, False
        

def get_dem_params(key='', file='./all_targets.pickle', all=True, plot=False, time_weighted=False, seconds_per=5, return_loc=False,
                regions_return=False, paramssaved=True, fromhome=False, keydict={},
                   doparam='above10s', namesearchstring=''):

    """
    set paramssaved=True to retrieve already-saved params in the dem result files 
        (viz.get_DEM_params() does this with the save_params_file attribute set True)
    """
    
    import glob

    if all:
        if fromhome:
            if not keydict:
                print('To get the result files from their homes, you need to set keydict equal to the dictionary containing key info.')
                return
                
            all_res_files=[]
            allkeys = keydict.keys()
            for kk in allkeys:
                res_files, tworegion = get_key_resultfiles(kk, file, withparams=paramssaved, fromhome=fromhome, keydict=keydict, 
                                                           namesearchstring=kk+'_'+namesearchstring)
                if tworegion:
                    all_res_files.extend(res_files[0])
                    all_res_files.extend(res_files[1])
                else:
                    all_res_files.extend(res_files)
        else:
                      
            if paramssaved:
                all_res_files = glob.glob('./compact_results/*withparams.pickle')
            else:
                rfs = glob.glob('./compact_results/*.pickle')
                all_res_files = [f for f in rfs if 'withparams' not in f]
            
    elif key:
        res_files, tworegion = get_key_resultfiles(key, file, withparams=paramssaved, fromhome=fromhome, keydict=keydict, 
                                                           namesearchstring=namesearchstring)

        if tworegion:
            res=[]
            if res_files[0]:
                res0 = extract_and_plot_param_histograms(res_files[0], key=key, doparam=doparam,
                                             plot=plot, time_weighted=time_weighted, seconds_per=seconds_per, plotlabel='Region 0')
                res.append(res0)
                
            if res_files[1]:
                res1 = extract_and_plot_param_histograms(res_files[1], key=key, doparam=doparam,
                                             plot=plot, time_weighted=time_weighted, seconds_per=seconds_per, plotlabel='Region 1')  
                res.append(res1)

            if regions_return:
                return res 

            all_res_files = res_files[0] + res_files[1]

        else:
            all_res_files = res_files
                
    else:
        print('Either set all=True, or select a key!')
        return

    #print(all_res_files)

    res = extract_and_plot_param_histograms(all_res_files, key=key, doparam=doparam, 
                                            plot=plot, time_weighted=time_weighted, seconds_per=seconds_per)

    if return_loc and key:
        with open(file, 'rb') as f:
            data = pickle.load(f)

        ARDict = data[key]
        loc = ARDict['loc']
        
        return res, loc
    else:
        return res



def extract_and_plot_param_histograms(res_files, key='', doparam='above10s',
                                      plot=False, time_weighted=False, seconds_per=5, 
                                      accthreshold=95,
                                      plotlabel='', show=False):


    import pandas as pd
    import astropy.time

    doishi=False
    dopeak=False
    dolow=False
    dohi=False

    flare_res = get_saved_flares(flarepath='./reference_files/', 
                                 add_stdv_flares=True, add_manual_flares=True)
    early_starts = flare_res[0]
    late_stops = flare_res[1]

    if doparam=='above10s':
        savekey = 'above_10MK'
        lowhigh=True
        colors1 = ['skyblue', 'aliceblue', 'lightcyan']
        colors2 = ['purple', 'lavender', 'thistle']
        thexlim = [1e17,1e25]
        doishi=True
        thexlabel='EM Integrated >10 MK'
        thexscale='log'

        thebins = np.geomspace(1e17, 1e25, 55)

    if doparam=='above5s':
        savekey = 'above_5MK'
        lowhigh=True
        colors1 = ['palevioletred', 'pink', 'lavenderblush']
        colors2 = ['purple', 'lavender', 'thistle']
        thexlim = [1e18,1e26]
        thexlabel='EM Integrated >5 MK'
        thexscale='log'

        thebins = np.geomspace(1e18, 1e26, 55)

    if doparam=='max_temp':
        savekey = 'max_temp'
        lowhigh=False
        colors1 = ['darkred', 'pink', 'lavenderblush']
        colors2 = ['purple', 'lavender', 'thistle']
        thexlim = [6.1,6.55]
        dopeak=True
        thexlabel='Peak Temperature'
        thexscale='linear'

        thebins = np.arange(6.1, 6.55, 0.01)

    if doparam=='low_powers':
        savekey = 'powerlaws'
        lowhigh=False
        colors1 = ['orange', 'khaki', 'lemonchiffon']
        colors2 = ['purple', 'lavender', 'thistle'] 
        thexlim=[1,4]
        dolow=True
        thexlabel='Lower Power Law'
        thexscale='linear'

        thebins = np.arange(1, 4, 0.05)

    if doparam=='hi_powers':
        savekey = 'powerlaws'
        lowhigh=False
        colors1 = ['red', 'khaki', 'lemonchiffon']
        colors2 = ['purple', 'lavender', 'thistle'] 
        thexlim=[-15,0]
        dohi=True
        thexlabel='Upper Power Law'
        thexscale='linear'

        thebins = np.arange(-15, 0, 0.3)

    
    
    all_params_flares = []
    all_params_non = []
    all_params=[]
    if lowhigh:
        all_paramsl_flares = []
        all_paramsl_non = []
        all_paramsl=[]
        all_paramsh_flares = []
        all_paramsh_non = []
        all_paramsh=[]   
        
    times = []
    durations =  []
    time1 = 0

    EMTall = []
    EMTflare = []
    EMTnon = []
    
    for f in res_files:
        #print(f)
        data, timestring, time = viz.load_DEM(f)
        #print(time[0])
        
        if savekey in data.keys():
            #print('found saved')
            param = data[savekey]
            chanax = data['chanax']
            EMT_thresh = data['EMT_thresh_5']
            dn_in = data['dn_in']

            
        else:
            #print('getting params')
            res = viz.get_DEM_params(f,save_params_file=True)
            if not res:
                print('couldnt get dem params.')
                continue
    
            #print(res)
        
            m1, max1, above5_, above7_, above10_, \
                above_peak, below_peak, above_635, below_635, \
                   chanax, dn_in, edn_in, \
                        powerlaws, EMT_all, EMT_thresh = res

            if doparam=='above10s':
                param = above10_

            if doparam=='above5s':
                param = above5_
                
            if doparam=='max_temp':
                param = m1
                #print(param)
                
            if doparam=='low_powers' or doparam=='hi_powers':
                param = powerlaws


        if doparam=='low_powers':
            param = param[0][0]
        if doparam=='hi_powers':
            param = param[1][0]



        #SHOULD IT BE HERE???
        #savefile=f.split('.p')[-2]+'_withparams.pickle'
        #data, timestring, time = viz.load_DEM(savefile)
        res = check_avg_rej(time, data['nustar_datapath'], threshold=accthreshold)
        if not res[1]:
            print('For time, ', time[0].strftime('%D %H-%M-%S'), '-', 
                  time[1].strftime('%D %H-%M-%S'), ' mean accepted events, ', 
                  res[0], ' below threshold, ', accthreshold)
            continue        
            
        times.append(time)
        dur = (time[1]-time[0]).to(u.s)
        if time1==0:
            time1=time[0]
        durations.append(dur)
        time_mult = round(dur.value/seconds_per)

        

        flare = check_for_flare(time, early_starts, late_stops)

        #if above10_[0] > 1e24:
        #    print(f, above10_)

        if len(dn_in) <= 6:
            print('Suspect missing nustar:')
            print(chanax)
            print(f)
            continue


    
        if flare:
            if time_weighted:

                if lowhigh:
                    all_params_flares.extend([param[0] for i in range(0,time_mult)])
                    all_paramsl_flares.extend([param[1] for i in range(0,time_mult)])
                    all_paramsh_flares.extend([param[2] for i in range(0,time_mult)])
                else:
                    all_params_flares.extend([param for i in range(0,time_mult)])
                EMTflare.extend([EMT_thresh for i in range(0,time_mult)])
            else:
                if lowhigh:
                    all_params_flares.append(param[0])
                    all_paramsl_flares.append(param[1])
                    all_paramsh_flares.append(param[2])
                else:
                    all_params_flares.append(param)
                EMTflare.append(EMT_thresh)
                
        else:
            if time_weighted:
                #print(time_mult)
                if lowhigh:
                    all_params_non.extend([param[0] for i in range(0,time_mult)])
                    all_paramsl_non.extend([param[1] for i in range(0,time_mult)])
                    all_paramsh_non.extend([param[2] for i in range(0,time_mult)])
                else:
                    all_params_non.extend([param for i in range(0,time_mult)])
                EMTnon.extend([EMT_thresh for i in range(0,time_mult)])
            else:  
                if lowhigh:
                    all_params_non.append(param[0])
                    all_paramsl_non.append(param[1])
                    all_paramsh_non.append(param[2])
                else:
                    all_params_non.append(param)
                EMTnon.append(EMT_thresh)
        
        if time_weighted:
            if lowhigh:
                all_params.extend([param[0] for i in range(0,time_mult)])
                all_paramsl.extend([param[1] for i in range(0,time_mult)])
                all_paramsh.extend([param[2] for i in range(0,time_mult)])
            else:
                all_params.extend([param for i in range(0,time_mult)])
            EMTall.extend([EMT_thresh for i in range(0,time_mult)])
        else:
            if lowhigh:
                all_params.append(param[0])
                all_paramsl.append(param[1])
                all_paramsh.append(param[2])
            else:
                all_params.append(param)
            EMTall.append(EMT_thresh)

    # #print(durations)
    # for d in durations:
    #     print(round(d.value/5))


    #print(all_params_flares)
    #print(all_params_non)
    #print(thebins)


    if plot:

        area_i = 100**2
        area_m = np.pi*150**2
        #print(area_i, area_m)
        factor = area_m/area_i
        #factor = 1
        
        from matplotlib import pyplot as plt 

        fig, axes = plt.subplots(2, 1, figsize=(15,6), tight_layout = {'pad': 1})

        ax=axes[0]
        if lowhigh:
            ax.hist(np.array(all_paramsl_non)*factor, bins=thebins, color=colors1[1], edgecolor='black', 
                    label='Uncertainty Lowbound', alpha=0.9, linestyle='dotted')
            ax.hist(np.array(all_paramsh_non)*factor, bins=thebins, color=colors1[2], edgecolor='black', 
                    label='Uncertainty Highbound', alpha=0.9, linestyle='dashed')
            
        ax.hist(np.array(all_params_non)*factor, bins=thebins, color=colors1[0], edgecolor='black', label='Non-flare time bins', alpha=0.9)

        if key:
            ax.set_title(key+' '+plotlabel)
        else:
            ax.set_title('All')

        ax=axes[1]
        if lowhigh:
            ax.hist(np.array(all_paramsl_flares)*factor, bins=thebins, color=colors2[1], edgecolor='black', 
                    label='Uncertainty Lowbound', alpha=0.9, linestyle='dotted')
            ax.hist(np.array(all_paramsh_flares)*factor, bins=thebins, color=colors2[2], edgecolor='black', 
                    label='Uncertainty Highbound', alpha=0.9, linestyle='dashed')
            
        ax.hist(np.array(all_params_flares)*factor, bins=thebins, color=colors2[0], edgecolor='black', label="Bins during flares")


        for ax in axes:
            ax.set_xscale(thexscale)
            ax.set_xlim(thexlim)
            ax.set_xlabel(thexlabel)
            
            if time_weighted:
                 ax.set_ylabel('# of '+str(seconds_per)+'s intervals')
            if doishi:
                ax.axvline(1.8e22, color='Red')
                ax.axvline(1.5e23, color='Red')
                ax.axvspan(1.8e22, 1.5e23, alpha=0.3, color='Red', label='Ishikawa (2017) 95% Interval')
                
            if dopeak:
                ax.axvline(6.505, color='Red')
                ax.axvline(6.602, color='Red')
                ax.axvspan(6.505, 6.602, alpha=0.3, color='Red', label='Warren (2012) Peaks')
    
                ax.axvline(6.301, color='cyan')
                ax.axvline(6.415, color='cyan')
                ax.axvspan(6.301, 6.415, alpha=0.2, color='cyan', label='Duncan (2024) Peaks')
                
            if dolow:
                ax.axvline(2, color='Red')
                ax.axvline(5, color='Red')
                ax.axvspan(2, 5, alpha=0.3, color='Red', label='Bradshaw (2012) Indices')

                ax.axvline(1.9, color='cyan')
                ax.axvline(2.3, color='cyan')
                ax.axvspan(1.9, 2.3, alpha=0.2, color='cyan', label='Duncan (2024) Indices')
                
            if dohi:
                ax.axvline(-10, color='Red')
                ax.axvline(-7, color='Red')
                ax.axvspan(-10, -7, alpha=0.3, color='Red', label='Barnes (2016)a Indices')

                ax.axvline(-12, color='cyan')
                ax.axvline(-4, color='cyan')
                ax.axvspan(-12, -4, alpha=0.2, color='cyan', label='Duncan (2024) Indices')
            
            if not key:
                if doparam=='above10s':
                    pass
                    #ax.set_ylim([0,4000])
                    #print('')
                if doparam=='above5s':
                    ax.set_ylim([0,4000])
                if doparam=='max_temp':
                    ax.set_ylim([0,7000])
                if doparam=='low_powers':
                    ax.set_ylim([0,3000])
                if doparam=='hi_powers':
                    ax.set_ylim([0,4000])

            ax.legend()

        
        if time_weighted:
            if key:
                plt.savefig('figures_etc/'+key+'_'+str(seconds_per)+'s_bin_time_weighted_'+doparam+'plot_'+plotlabel+'.png')
            else:
                plt.savefig('figures_etc/all_'+str(seconds_per)+'s_bin_time_weighted_'+doparam+'plot_'+plotlabel+'.png')
        else:
            if key:
                plt.savefig('figures_etc/'+key+'_'+doparam+'_plot_'+plotlabel+'.png')
            else:
                plt.savefig('figures_etc/all_'+doparam+'_plot_'+plotlabel+'.png')

        if not show:
            plt.close()


    EMT = [EMTall, EMTflare, EMTnon]

    values = [all_params, all_params_flares, all_params_non]
    averages = [np.mean(all_params), np.mean(all_params_flares), np.mean(all_params_non)]
    
    if lowhigh:
        averagesl = [np.mean(all_paramsl), np.mean(all_paramsl_flares), np.mean(all_paramsl_non)]
        averagesh = [np.mean(all_paramsh), np.mean(all_paramsh_flares), np.mean(all_paramsh_non)]
        return values, np.array(averages), np.array(averagesl), np.array(averagesh), time1, EMT

    else:
        return values, np.array(averages), time1, EMT

    
def check_for_flare(time, starts, stops):

    flare=False
    #If there is a flare that starts before time0 and stops after time0
    b4 = [s < time[0] for s in starts]
    ea = [s > time[0] for s in stops]
    es = np.where(np.logical_and(b4, ea))

    if es[0].size > 0:
        flare=True

    #If there is a flare that starts between time0 and time1
    b4 = [s > time[0] for s in starts]
    ea = [s < time[1] for s in starts]
    es = np.where(np.logical_and(b4, ea))

    if es[0].size > 0:
        flare=True

    return flare


def get_saved_flares(flarepath='./', add_stdv_flares=False, add_manual_flares=True,
                      specific_key_file=[]):

    
    if specific_key_file:
        with open(specific_key_file, 'rb') as f:
            data = pickle.load(f)
        early_starts, late_stops = data['early_starts'], data['late_stops']
        
        return [early_starts, late_stops]
        
    import pandas as pd
    import astropy.time as time

    risefactor=0.5
    fallfactor=2
    
    df = pd.read_csv(flarepath+'fpmA.csv')
    starts = df['flare_start'].values
    stops = df['flare_end'].values
    num=len(starts)
    
    from astropy import units as u
    durs = [(time.Time(stops[i])-time.Time(starts[i])).to(u.s) for i in range(0,num)]
    early_starts = [(time.Time(starts[i])-risefactor*durs[i]).datetime for i in range(0,num)]
    late_stops = [(time.Time(stops[i])+fallfactor*durs[i]).datetime for i in range(0,num)]

    #early_starts = [(time.Time(s)-2*u.min) for s in starts]
    #late_stops = [(time.Time(s)+4*u.min) for s in stops]

    if add_stdv_flares:
        with open(flarepath+'stdv_flares.pickle', 'rb') as f:
            data = pickle.load(f)

        flarewindows = np.array(data['stdv_flares'])
        num=len(flarewindows)
        durs = [time.Time(flarewindows[i,1])-time.Time(flarewindows[i,0]) for i in range(0, num)]
        early_starts.extend([(time.Time(flarewindows[i,0])-risefactor*durs[i]).datetime for i in range(0, num)])
        late_stops.extend([(time.Time(flarewindows[i,1])+fallfactor*durs[i]).datetime for i in range(0, num)])
        #early_starts.extend([(time.Time(s)-2*u.min).datetime for s in flarewindows[:,0]])
        #late_stops.extend([(time.Time(s)+4*u.min).datetime for s in flarewindows[:,1]]) 


    if add_manual_flares:
        with open(flarepath+'manual_flares.pickle', 'rb') as f:
            data = pickle.load(f)

        #Not adding a rise/fall factor bc manual inspection looked for whole flare times.
        flarewindows = np.array(data['manual_flares'])
        num=len(flarewindows)
        early_starts.extend([flarewindows[i,0].datetime for i in range(0, num)])
        late_stops.extend([flarewindows[i,1].datetime for i in range(0, num)])   

    return [early_starts, late_stops]





def make_orbit_plots(working_dir, key, minT=5.6, maxT=7.2, show=False):

    
    all_time_intervals, all_time_intervals_list = tis.find_all_intervals(working_dir, shush=True, 
                                                                        missing_last=False)
    print(len(all_time_intervals))
    
    
    print('Number of orbits:', len(all_time_intervals))
    for k in range(0, len(all_time_intervals)):
        time_intervals = all_time_intervals[k]
    

        vals = viz.get_DEM_timeseries(time_intervals, working_dir, minT, maxT, key)    

        if len(vals['result_time_intervals'])==0:
            print('No sucessful DEM intervals.')
            continue
        
        
        dn_inss=vals['dn_ins']
        labels= ['NuSTAR 6-10 keV Emission', 'NuSTAR 3.5-6 keV Emission', 'NuSTAR 2.5-3.5 keV Emission']
        backcolors=['turquoise', 'lightcyan']
        color='darkcyan'
        for i in range(1,4):
            ind = i*-1
            nuvals = [din[ind] for din in dn_inss]    
            viz.pretty_orbit_timeseries(vals['result_time_intervals'], nuvals, 'Cts/s', labels[i-1],
                                color, backcolors, working_dir=working_dir, plot_flares=True, show=show)
        
        
        peaks=vals['peaks']
        peaksmk = [10**m1/1e6 for m1 in peaks]    
        
        backcolors=['pink', 'lavenderblush']
        color='Red'
            
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], peaksmk, 'DEM Peak Temperature (MK)', 'DEM Peak Temperature',
                                color, backcolors, working_dir=working_dir, show=show)
        
        
        backcolors=['powderblue', 'aliceblue']
        color='Blue'
        
        above10s=np.array(vals['above10s'])
        above10s_=above10s[:,0]
        
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], above10s_, 'EM (cm^-5)', 'Total EM >10 MK',
                                color, backcolors, error=True, quantity_low=above10s[:,1], quantity_high=above10s[:,2], 
                                ylog=True, comparisonbar=True, comp_band=[1.8e22, 1.5e23, 'Ishikawa (2017) 95%'],
                                    working_dir=working_dir, plot_flares=True, show=show)
        
        backcolors=['powderblue', 'aliceblue']
        color='Green'
        
        above7s=np.array(vals['above7s'])
        above7s_=above7s[:,0]
        
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], above7s_, 'EM (cm^-5)', 'Total EM >7 MK',
                                color, backcolors, error=True, quantity_low=above7s[:,1], quantity_high=above7s[:,2], 
                                ylog=True, working_dir=working_dir, show=show)
        
        backcolors=['powderblue', 'aliceblue']
        color='Purple'
        
        above5s=np.array(vals['above5s'])
        above5s_=above5s[:,0]
        
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], above5s_, 'EM (cm^-5)', 'Total EM >5 MK',
                                color, backcolors, error=True, quantity_low=above5s[:,1], quantity_high=above5s[:,2], 
                                ylog=True, working_dir=working_dir, show=show)
        
        
        backcolors=['khaki', 'lemonchiffon']
        color='Orange'
        
        val=np.array(vals['low_powers'])
        
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], val, 'Index', 'Lower Power Law',
                                color, backcolors, error=False, working_dir=working_dir, show=show)
        
        
        backcolors=['khaki', 'lemonchiffon']
        color='Red'
        
        val=np.array(vals['hi_powers'])*-1
        
        viz.pretty_orbit_timeseries(vals['result_time_intervals'], val, 'Index', 'Upper Power Law',
                                color, backcolors, error=False, working_dir=working_dir, show=show)






def get_region_directories(key, targets_file='./all_targets.pickle'):

    import glob
    import os

    with open(targets_file, 'rb') as f:
        data = pickle.load(f)

    ARDict = data[key]
    working_dir = ARDict['working_dir']

    method = ARDict['method']
    #print(method)
    
    if method=='fit':
        directories=[working_dir]


    if method=='input':
        directories = [f+'/' for f in glob.glob(working_dir+'/region_*') if os.path.isdir(f)]
        directories.sort()

    if method=='double':
        gauss_stats = ARDict['gauss_stats']
        sep_axis = gauss_stats[0][0]
        if sep_axis=='EW':
            directions = ['east', 'west']
        elif sep_axis=='SN':
            directions = ['south', 'north']
            
        directories=[]    
        for r in directions:
            directories.append(working_dir+'/'+r+'/')

    return directories


def do_stdv_analysis(key, file, show=True):

    import nustar_utilities as nuutil
    import glob

    with open(file, 'rb') as f:
        data = pickle.load(f)

    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    method = ARDict['method']
    directories = get_region_directories(key, targets_file=file)

    all_newwindows=[]

    for dd in directories:
    
        all_time_intervals, all_time_intervals_list = tis.find_all_intervals(dd, shush=True, 
                                                                            missing_last=False)
    
        for ind in range(0, len(all_time_intervals)):
            datapath=id_dirs[ind]
            obsid=obsids[ind]
            print(datapath)

            evt_data, hdr = nu.return_submap(datapath=datapath, fpm='A', return_evt_hdr=True)
            time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
            timerange = [time0.tt.datetime, time1.tt.datetime]
            from datetime import timezone
            
            #Comment second line if you're not using this same example nustar orbit
            #Edit it to include only the desired time interval (default- all times in file) once you've run this once
            #timerange=[]
            #timerange=[datetime.datetime(2018, 5, 29, 16, 00), datetime.datetime(2018, 5, 29, 16, 55)]
            timerange = [t.replace(tzinfo=timezone.utc) for t in timerange]


            evtA = glob.glob(datapath+'/event_cl/*A06_cl.evt')
            evtB = glob.glob(datapath+'/event_cl/*B06_cl.evt')
            hkA  = glob.glob(datapath+'/hk/*A_fpm.hk')
            hkB  = glob.glob(datapath+'/hk/*B_fpm.hk')

            import lightcurves as lc
            
            lc.prepare_nustar_lightcurves(evtA, evtB, hkA, hkB, timebin=15, erange=[2.,4.], 
                                          livetime_corr=False, save_dir=dd, event_stats=True)
            lc.prepare_nustar_lightcurves(evtA, evtB, hkA, hkB, timebin=15, erange=[4.,6.], 
                                          livetime_corr=False, save_dir=dd, event_stats=True)
            lc.prepare_nustar_lightcurves(evtA, evtB, hkA, hkB, timebin=15, erange=[6.,10.], 
                                          livetime_corr=False, save_dir=dd, event_stats=True)
            
            newwindows = lc.plot_with_stdv(aia_inputs=[], fexviii=False, nustar_inputs=[[2.,4.],[4.,6.],[6.,10.]], 
                   goes_inputs=[], smooth=24, 
                   plot_each=True, plot_logic=True, remove_fexviii_max=False, 
                   analyze_transients=True, transient_number=3,
                   timerange=timerange,
                  excluded_range=[], save_dir=dd, savestring=key+'_'+obsid, show=show)
            

            all_newwindows.extend(newwindows)


    return all_newwindows





def make_summary_lcs(key, file, flarepath='./reference_files/', specific_region_flares=False,
                     #method='input', 
                     show=True, goes=True,
                    accthreshold=95, pre_dem_nustar_only=False,
                    use_inputs_only=False):

    
    from matplotlib import pyplot as plt
    import copy
    import matplotlib.dates as mdates
    import nustar_utilities as nuutil
    import lightcurves as lc
    import glob

    minT=5.6
    maxT=7.2
    nustar_acc_color='xkcd:apple'
    nustar_cts_color=['xkcd:sky blue', 'xkcd:cerulean', 'xkcd:periwinkle', 'xkcd:cyan']

    with open(file, 'rb') as f:
        data = pickle.load(f)

    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    directories = get_region_directories(key, targets_file=file)

    for dd in directories:
    
        all_time_intervals, all_time_intervals_list = tis.find_all_intervals(dd, shush=True, 
                                                                            missing_last=False)
    
        for ind in range(0, len(all_time_intervals)):
            
            datapath=id_dirs[ind]
        
            time_intervals = all_time_intervals[ind]
            if not pre_dem_nustar_only:
                vals = viz.get_DEM_timeseries(time_intervals, dd, minT, maxT, key,
                                               inputs_only=use_inputs_only)   
                time_intervals = vals['result_time_intervals']
            
            if not time_intervals:
                print('This key/region/orbit had no sucessful DEMs.')
                print(key)
                print(dd)
                print(ind)
                #print('')
                return
            
            evt_data, hdr, obsid = nu.return_submap(datapath=datapath, fpm='A', return_evt_hdr=True, return_obsid=True)
            time0, time1 = [nuutil.convert_nustar_time(hdr['TSTART']), nuutil.convert_nustar_time(hdr['TSTOP'])]
            timerange = [time0.datetime, time1.datetime]
            print(timerange[0].strftime('%H-%M-%S'), timerange[1].strftime('%H-%M-%S'))
            
            
            from datetime import timezone
            newtimerange = [t.replace(tzinfo=timezone.utc) for t in timerange]
            
            axiscolor='black'
            
            labelfontsize=10
            tickfontsize=17
            
            
            fig, axes = plt.subplots(6, 1, figsize=(15, 15), sharex=True)
            plt.subplots_adjust(hspace=0)
            
            
            if not pre_dem_nustar_only:
                # AIA PLOTS + NuSTAR DEM INPUT PLOTS + >10 MK EM plots #======================================================================
                
                lw=2
                times = [t[0].datetime for t in time_intervals]
                midtimes=[(t[0]+(t[1]-t[0]).to(u.s)/2).datetime for t in time_intervals]
                
                
                starttime = (time_intervals[0][0]-120*u.s).datetime
                stoptime = (time_intervals[-1][1]+120*u.s).datetime
                
                times_ = copy.deepcopy(times)
                times_.append(time_intervals[-1][1].datetime)
                
                
                dn_ins = vals['dn_ins']
                chanaxs = vals['chanaxs']
                
                allcolors = lc.make_colors(10)
                
                if not use_inputs_only:

                    ax=axes[0]
                    normin=1
                    for i in range(0, 6):  
                        aiavals = [din[i] for din in dn_ins] 
                        label = chanaxs[0][i]
                        color = allcolors[i]
                        normvals = np.array(aiavals)/np.max(aiavals)
                        if np.min(normvals) < normin:
                            normin=np.min(normvals)
                        ax.stairs(normvals, times_, linewidth=lw, color=color, 
                                  label=label,
                                   baseline=None)
                    ax.set_ylim([normin*0.95, 1.01])
                
                ax=axes[1]
                normin=1
                for i in range(0, 3): 
                    ii = (i+1)*-1
                    aiavals = [din[ii] for din in dn_ins] 
                    label = chanaxs[0][ii]
                    color = allcolors[ii]
                    normvals = np.array(aiavals)/np.max(aiavals)
                    if np.min(normvals) < normin:
                        normin=np.min(normvals)
                    ax.stairs(normvals, times_, linewidth=lw, color=color, 
                              label=label,
                             baseline=None)
                
                ax.set_ylim([normin*0.95, 1.01])
                
                ax=axes[5]
                if not use_inputs_only:
                    ax=axes[5]
                    above10s=np.array(vals['above10s'])
                    above10s_=above10s[:,0]
                    label='Total EM >10 MK'
                    ylabel='EM (cm^-5)'
                    color='Blue'
                    ax.stairs(above10s_, times_, linewidth=lw, color=color, 
                              label=label,
                             baseline=None)
                
                    quantity_low=above10s[:,1]
                    quantity_high=above10s[:,2]
                
                    error=True
                    if error:
                        lp = np.hstack([quantity_low, quantity_low[-1]])
                        hp = np.hstack([quantity_high, quantity_high[-1]])
        
                        fill = ax.fill_between(times_, lp, hp, step="post", 
                                         color=color, alpha=0.1) 
                
                    ax.set_ylabel(ylabel)
                    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                    ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))
                    ax.set_yscale('log')
                
                comp_band=[1.8e22, 1.5e23, 'Ishikawa (2017) 95%']
                comparisonbar=True
                if comparisonbar:
                    ax.axhspan(comp_band[0], comp_band[1], color='brown', alpha=0.5, label=comp_band[2])   
                    ax.legend()
                
                
            
            #==================================================================================================================================
            # Mark times where mean accepted events level over course of time interval is below threshold

            for tr in time_intervals:
                meanevsum, checkacc = check_avg_rej(tr, datapath, threshold=accthreshold)
                if not checkacc:
                    midtime = tr[0] + (tr[1]-tr[0]).to(u.s).value/2*u.s
                    for ax in axes:
                        ax.axvline(midtime.datetime, color='aqua', lw=10, zorder=-1)



            #=================================================================================================================================
            
            #NUSTAR - ALL-FOV, FINE TIME RESOLUTION ==========================================================================================
            
            tr=newtimerange
                
            evtA = glob.glob(datapath+'/event_cl/*A06_cl.evt')
            evtB = glob.glob(datapath+'/event_cl/*B06_cl.evt')
            hkA  = glob.glob(datapath+'/hk/*A_fpm.hk')
            hkB  = glob.glob(datapath+'/hk/*B_fpm.hk')
            
            #Load in the evt file (has the list of photons)
            evtdataA, hdrA = lc.load_nufiles(evtA[0])
            # Load in the hk file (has the livetime info)
            lvdataA, lvhdrA = lc.load_nufiles(hkA[0])
            evtdataB, hdrB = lc.load_nufiles(evtB[0])
            lvdataB, lvhdrB = lc.load_nufiles(hkB[0])
    
            
            eranges=[[2.5,3.5],[3.5,6.], [6.,10.], [10.,15.]]
            labels = ['NuSTAR 2.5-3.5 keV', 'NuSTAR 3.5-6 keV', 'NuSTAR 6-10 keV', 'NuSTAR 10-15 keV']
            eind=0
            #     acc=0
            for erange in eranges:
                kevA = evtdataA['PI']*0.04+1.6
                erange_evtdataA = evtdataA[np.where(np.logical_and(kevA > erange[0],kevA < erange[1]))]
                kevB = evtdataB['PI']*0.04+1.6
                erange_evtdataB = evtdataB[np.where(np.logical_and(kevB > erange[0],kevB < erange[1]))]
            
            
                res = lc.get_a_nustar_lightcurve(erange_evtdataA, hdrA, lvdataA, lvhdrA, timebin=10, livetime_corr=True, event_stats=True)
                times_converted, countrate, lvt, counts, acc_sample, rej_sample, all_sample = res
              
                keepinds = np.nonzero(np.logical_and(times_converted > tr[0], times_converted < tr[1]))
                times_converted=times_converted[keepinds]
                countrateA=countrate[keepinds]   

                acc_sample=acc_sample[keepinds]
                rej_sample=rej_sample[keepinds]
                evsumA=acc_sample/(acc_sample+rej_sample)*100

                    
                res = lc.get_a_nustar_lightcurve(erange_evtdataB, hdrB, lvdataB, lvhdrB, timebin=10, livetime_corr=True, event_stats=True)
                times_converted, countrate, lvt, counts, acc_sample, rej_sample, all_sample = res
    
            
                keepinds = np.nonzero(np.logical_and(times_converted > tr[0], times_converted < tr[1]))
                times_converted=times_converted[keepinds]
                countrateB=countrate[keepinds]

                acc_sample=acc_sample[keepinds]
                rej_sample=rej_sample[keepinds]
        
                evsumB=acc_sample/(acc_sample+rej_sample)*100


                evsum = (np.array(evsumA)+np.array(evsumB))/2.

                
                axes[3].plot(times_converted, evsum, label='Accepted Event %', color=nustar_acc_color)
    
                totalrate = countrateA+countrateB
                totalrate_vals = totalrate[np.isfinite(totalrate)]
                maximum = np.nanmax(totalrate_vals) 
                axes[2].plot(times_converted, totalrate/maximum, label=labels[eind], color=nustar_cts_color[eind]) 

            
                eind+=1
                
            axes[2].set_ylim([0,1.01])
            
            #==================================================================================================================================
            
            #Add Goes Curve (both days) #=========================================================================================================
            
            if goes:
                tr=timerange
    
                goes_number = ARDict['goes_satellite']
                
                #print(tr)
                lc.get_goes(tr, satellite=goes_number)
                instrument='GOES'
                data = lc.load_lightcurves(instrument)
                
                ylabel = data['GOES flux label']
                goestimes = data['GOES Times']
                xrsbcounts = data['XRSB counts']
                xrsblabel = data['XRSB label']
                gts = [t.datetime for t in goestimes]
                
                axg = axes[4] #.twinx()
                #xrs=[x.value for x in all_xrsb]
                axg.plot(gts, xrsbcounts, label='GOES '+xrsblabel, color='red')
                axg.set_yscale('log')
                #axg.set_ylim([1e-7,1e-6])
                #axg.yaxis.set_ticks([1e-6, 1e-5, 1e-4], labels=['C','M','X'], color='red', fontsize=labelfontsize)
                axg.grid(True, which='major', axis='y', linestyle=':', color='red')
                
                #axg.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                #axg.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))
            
            #==================================================================================================================================
            
            # #Plot Flares from Reed's list #====================================================================================================
            # plot_flares=True
            # if plot_flares:
            #     import pandas as pd
            #     import astropy.time
            #     df = pd.read_csv('fpmA.csv')
            #     starts = df['flare_start'].values
            #     stops = df['flare_end'].values
                
            
            #     early_starts = [(astropy.time.Time(s)-2*u.min).datetime for s in starts]
            #     late_stops = [(astropy.time.Time(s)+4*u.min).datetime for s in stops]  
            
            #     for j in range(0, len(early_starts)):
            #         for ax in axes:
            #             ax.axvspan(early_starts[j], late_stops[j], alpha=.25, color='purple')
            #             ax.set_xlim(timerange[0], timerange[1])


            #Plot "Flares" from Reed's other list #========================================================================================
            # plot_flares=True
            # if plot_flares:
            #     import pandas as pd
            #     import astropy.time
            #     df = pd.read_csv('fpmA_notflares.csv')
            #     starts = df['start'].values
            #     stops = df['end'].values
                
            
            #     early_starts = [(astropy.time.Time(s)-1*u.min).datetime for s in starts]
            #     late_stops = [(astropy.time.Time(s)+1*u.min).datetime for s in stops]  
            
            #     for j in range(0, len(early_starts)):
            #         for ax in axes:
            #             ax.axvspan(early_starts[j], late_stops[j], alpha=.5, color='cyan')
            #             ax.set_xlim(timerange[0], timerange[1])
            
            

            # #Plot Flares stdv analysis list #====================================================================================================
            # plot_stdv_flares=True
            # if plot_stdv_flares:
            #     #print('got here')

            #     with open('stdv_flares.pickle', 'rb') as f:
            #         data = pickle.load(f)

            #     flarewindows = np.array(data['stdv_flares'])
                
            #     early_starts = [(astropy.time.Time(s)-2*u.min).datetime for s in flarewindows[:,0]]
            #     late_stops = [(astropy.time.Time(s)+4*u.min).datetime for s in flarewindows[:,1]]  
            
            #     for j in range(0, len(early_starts)):
            #         for ax in axes:
            #             ax.axvspan(early_starts[j], late_stops[j], alpha=.25, color='green')
            #             ax.set_xlim(timerange[0], timerange[1])
                        
            
            plot_flares=True
            if not specific_region_flares:
                plot_stdv_flares=True
                specific_key_file=[]
            else:
                plot_stdv_flares=False
                specific_key_file=ARDict['flarefile']
            
            if plot_flares:
                flare_res = get_saved_flares(flarepath=flarepath, specific_key_file=specific_key_file)
                early_starts = flare_res[0]
                late_stops = flare_res[1]

                for j in range(0, len(early_starts)):
                    for ax in axes:
                        ax.axvspan(early_starts[j], late_stops[j], alpha=.25, color='blue')
                        ax.set_xlim(timerange[0], timerange[1])

                if plot_stdv_flares:
                    flare_res = get_saved_flares(flarepath=flarepath, 
                                                 add_stdv_flares=True, add_manual_flares=True)
                    early_starts = flare_res[0]
                    late_stops = flare_res[1]
                    for j in range(0, len(early_starts)):
                        for ax in axes:
                            ax.axvspan(early_starts[j], late_stops[j], alpha=.25, color='red')
                            ax.set_xlim(timerange[0], timerange[1])
                    

            
            
            #==================================================================================================================================


            #Plot SAAs and Pointing Shifts #====================================================================================================
            plot_badtimes=False
            if plot_badtimes:
                import use_correlator as uc
                import astropy.time
                t_corr=2
                #minimum length of suborbit
                min_t = 30
                #minimum number of steps needed in a sufficiently-sized suborbit, based on the above. 
                min_step = int(np.ceil(min_t/t_corr))
                res = uc.get_suborbits(id_dirs[ind], t_corr, min_step, plot=False)
                bad_grouptimes = res[3]
                
            
                for j in range(0, len(bad_grouptimes)):
                    for ax in axes:
                        ax.axvspan(bad_grouptimes[j][0], bad_grouptimes[j][1], alpha=.25, color='orange')
                        ax.set_xlim(timerange[0], timerange[1])
            
            #==================================================================================================================================
            
            
            leg0 = axes[0].legend()
            leg1 = axes[1].legend()#(loc='best')
            leg2 = axes[2].legend()#loc='best')
            leg3 = axes[3].legend()
            leg4 = axes[4].legend()
            leg5 = axes[5].legend()
            
            for ax in axes:
                ax.set_xlim(timerange[0], timerange[1])
            
            
            plt.savefig(dd+obsid+'_lc_summary_new.png', transparent=False)

            if not show:
                plt.close()
                
        #print('')


def check_avg_rej(tr, datapath, threshold=95):
    """
    For the input time range (which must be during the interval covered by the data
    at location datapath), return the mean % of accepted events (meanevsum), as well as a 
    boolean regarding whether it is above or below the chosen threshold. 
    
    """
    import glob
    import lightcurves as lc
    
        
    evtA = glob.glob(datapath+'/event_cl/*A06_cl.evt')
    evtB = glob.glob(datapath+'/event_cl/*B06_cl.evt')
    hkA  = glob.glob(datapath+'/hk/*A_fpm.hk')
    hkB  = glob.glob(datapath+'/hk/*B_fpm.hk')
    
    #Load in the evt file (has the list of photons)
    evtdataA, hdrA = lc.load_nufiles(evtA[0])
    # Load in the hk file (has the livetime info)
    lvdataA, lvhdrA = lc.load_nufiles(hkA[0])
    
    evtdataB, hdrB = lc.load_nufiles(evtB[0])
    lvdataB, lvhdrB = lc.load_nufiles(hkB[0])
    
    
    erange=[2.5,3.5] #,[3.5,6.], [6.,10.], [10.,15.]]

    kevA = evtdataA['PI']*0.04+1.6
    erange_evtdataA = evtdataA[np.where(np.logical_and(kevA > erange[0],kevA < erange[1]))]
    kevB = evtdataB['PI']*0.04+1.6
    erange_evtdataB = evtdataB[np.where(np.logical_and(kevB > erange[0],kevB < erange[1]))]

    res = lc.get_a_nustar_lightcurve(erange_evtdataA, hdrA, lvdataA, lvhdrA, timebin=10, livetime_corr=True, event_stats=True)
    times_converted, countrate, lvt, counts, acc_sample, rej_sample, all_sample = res
  
    keepinds = np.nonzero(np.logical_and(times_converted > tr[0], times_converted < tr[1]))

    acc_sample=acc_sample[keepinds]
    rej_sample=rej_sample[keepinds]
    evsumA=acc_sample/(acc_sample+rej_sample)*100


    #print(len(evsumA))
        
    res = lc.get_a_nustar_lightcurve(erange_evtdataB, hdrB, lvdataB, lvhdrB, timebin=10, livetime_corr=True, event_stats=True)
    times_converted, countrate, lvt, counts, acc_sample, rej_sample, all_sample = res


    #keepinds = np.nonzero(np.logical_and(times_converted > tr[0], times_converted < tr[1]))

    acc_sample=acc_sample[keepinds]
    rej_sample=rej_sample[keepinds]

    evsumB=acc_sample/(acc_sample+rej_sample)*100


    #print(len(evsumB))
    
    evsum = (np.array(evsumA)+np.array(evsumB))/2.
    
    meanevsum = np.mean(evsum)

    if meanevsum <= threshold:
        return (meanevsum, False)
    else:
        return (meanevsum, True)


def plot_with_discriminator(axes, resfiles, discriminator, options, label, paramx, paramy, optionnumbers, caselabels,
                  xindex='', yindex='', flare=False, startstop=[], update_flarecheck=False):

    """

    Takes in plot axes, a set of result files from a given key/region, and selections re what parameters we would like to plot.

    paramx, paramy - to be plotted on the horizontal and vertical axes, respectively. Should be keys in the result file dictionaries. 

    xindex, yindex - for result params that come in lists/arrays (say, [value, lowbound, highbound]), indicates which to use. If set to
                        something other than an integer, will not be used (assume single value). 

    label - uniquely identifies the key/region; used for plot labels.

    Also takes in DISCRIMINATOR – a value of the key/region with which we would like to separate out different cases and plot them
    on different axes. The number of axes should be equal to the number of possible values of discriminator. For recordkeeping, we 
    also need: 
    
    optionnumbers – array (of the same length as the number of possible values of the discriminator), containing the number of each which
                    has been found so far (we'll add to the number associated with the value of discriminator each time this is run). 

    caselabels - list (of the same length as the number of possible values of the discriminator) containting lists of the LABELS of each
                    key/region (we'll add the present label value to the list for the present value of the discriminator.)
    
    If startstop is not set, take all intervals.
    If startstop contains lists of paired start and stop times for flares, then:
        if Flare == True: take only times intersecting with a flare
        if Flare == False: take only times that do not intersect with a flare

    update_flarecheck - set True to reevaluate the resultfiles to check if they are for time intervals intersecting with a flare. This
                            should be done if the flarelist changes. 
                            
    """

    from astropy import units as u

    greenz = ['green', 'darkgreen', 'mediumseagreen', 'lightgreen', 'forestgreen', 'seagreen', 'darkseagreen', 'lime', 'springgreen', 
          'mediumaquamarine', 'aquamarine', 'lawngreen', 'green', 'darkgreen', 'mediumseagreen', 'lightgreen', 'forestgreen', 'seagreen', 
              'darkseagreen', 'lime', 'springgreen', 
          'mediumaquamarine', 'aquamarine', 'lawngreen']

    bluez = ['blue', 'royalblue', 'mediumblue', 'darkblue', 'navy', 'cornflowerblue', 'dodgerblue', 'deepskyblue', 'skyblue',
         'lightblue', 'lightskyblue', 'steelblue', 'lightsteelblue', 'powderblue', 'darkslateblue', 'slateblue', 'lightslategray',
         'slategray']


    reds = ['lightcoral', 'indianred', 'brown', 'red', 'salmon', 'tomato']

    purples = ['plum', 'violet', 'purple', 'darkviolet', 'darkorchid', 'thistle']

    colorz = [greenz, bluez, reds, purples]

    first = 0

    #If we're using a list of flaretimes to decide whether to plot contents of these resultfiles...
    if startstop:
        early_starts = startstop[0]
        late_stops = startstop[1]

    valuex=[]
    valuey=[]   
    #For each result...
    for r in resfiles:
        #load in file
        with open(r, 'rb') as f:
            data = pickle.load(f)

        #If we're checking for a flare...
        if startstop:
            time = data['time_interval']
            #If there's already a note about whether it's a flare in the result file (and we aren't re-checking)
            if 'flare?' in data.keys() and update_flarecheck == False:
                flarecheck = data['flare?'] 
            else:
                #Check if this time interval intersects with a flare
                flarecheck = check_for_flare(time, early_starts, late_stops)
                #If set, update the result dictionary to make a note of the result.
                if update_flarecheck:
                    data['flare?'] = flarecheck
                    with open(r, 'wb') as f:
                         # Pickle the 'data' dictionary using the highest protocol available.
                         pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 
                        
            #If this time interval does not satisfy our condition re flaring, don't plot.            
            if flarecheck != flare:
                #print(time[0].strftime('%H-%M-%S'))
                continue
                
        #Check if x,yindex set; in either case, fetch the values needed.  
        if type(xindex) == int:
            valuex.append(data[paramx][xindex])
        else:
            valuex.append(data[paramx])
            
        if type(yindex) == int:
            valuey.append(data[paramy][yindex])
        else:
            valuey.append(data[paramy])

    #Check to see if the discriminator is equal to any of the options.
    for o in range(0, len(options)):
        #If so, plot it in the corresponding place.
        if discriminator == options[o]:
            ax=axes[o]
            #If this is the first time we're plotting a result from this specific key/region, add the label to our list.
            if label not in caselabels[o]:
                optionnumbers[o]+=1
                caselabels[o].append(label)
                first = 1
            thecolor = colorz[o][int(optionnumbers[o])]
            #If this is the first time we're plotting a result from this specific key/region, plot with label
            if first:
                ax.scatter(valuex, valuey, color=thecolor, label=label)
                ax.set_title(options[o]+' regions, '+paramy+' vs. '+paramx)
                firstl=0
            else:
                ax.scatter(valuex, valuey, color=thecolor)             


    return optionnumbers, caselabels



def consistency_hist_plot(all_sumcons, all_sumcons_non, all_sumcons_flare, dir_, label_='', time_weighted=False, ylabel='# of intervals',
                  show=False):

    import matplotlib.pyplot as plt

    tickfont=12
    ls=2.5

    fig, ax = plt.subplots(1, 1, figsize=(15,4), tight_layout = {'pad': 1})

    ax.hist(all_sumcons, bins=[0.5,1.5,2.5,3.5, 4.5, 5.5], edgecolor='black', color='palevioletred')
    ax.set_ylabel(ylabel)   
    ax.set_xlim([0, 6.5])
    ax.set_xticks([1,2,3,4,5,6], ['5', '7', '10', '12.6', '15.8', '(>15.8)'], fontsize=tickfont)
    ax.set_xlabel('Upper Limit (MK)', fontsize=15)
    ax.set_title(label_+' Histogram, all-orbits: lowest-temperature upper bound where DEM is consistent with NuSTAR 6-10 keV emission',
                     fontsize=15)
    plt.savefig(dir_+'/'+label_+'_all_orbit_consistency_histogram.png')
    if not show:
        plt.close()

    fig, axes = plt.subplots(2, 1, figsize=(15,8), tight_layout = {'pad': 1})

    ax=axes[0]
    ax.hist(all_sumcons_non, bins=[0.5,1.5,2.5,3.5, 4.5, 5.5], edgecolor='black', color='bisque', label='Non-flaring times')
    ax.set_title(label_+' Histogram, all-orbits: lowest-temperature upper bound where DEM is consistent with NuSTAR 6-10 keV emission',
                     fontsize=15)

    ax=axes[1]
    ax.hist(all_sumcons_flare, bins=[0.5,1.5,2.5,3.5, 4.5, 5.5], edgecolor='black', color='lavender', label='Flaring times')

    for ax in axes:
        if label_=='All_keys_':
            ax.set_ylim([0,6000])
        ax.set_ylabel(ylabel) 
        ax.set_xlim([0, 6.5])
        ax.set_xticks([1,2,3,4,5,6], ['5', '7', '10', '12.6', '15.8', '(>15.8)'], fontsize=tickfont)
        ax.set_xlabel('Upper Limit (MK)', fontsize=15)
        ax.legend(fontsize=tickfont)  
            

    plt.savefig(dir_+'/'+label_+'_flare_separated_all_orbit_consistency_histogram_rescale.png')
    if not show:
        plt.close()




    
def plot_temp_consistency(key='', file='all_targets.pickle', time_weighted=True, seconds_per=5, show=False, 
                         fetch_for_all=False, accthreshold=95, return_region_quiet_max=False):

    """
    Retrieves different DEM result files for runs with different temperature bounds, and then makes consistency plots:
        "Lowest-temperature upper bound where DEM is consistent with NuSTAR 6-10 keV emission"

    Set up to fetch result files from their time interval directories, and to work with all region selection methods.

    Keywords:

    fetch_for_all - set True to return all, flare, non arrays of values (time weighted or not) for incorporation
                    into a larger plot later (not made here). 

    key - date key, set to make plots for the regions associated with that key.
    time_weighted - set True to normalize histograms by time elapsed (vs. # of intervals)
    seconds_per - time interval used in time weighting (if time_weighted == True).
    
    
    """

    import time_interval_selection as tis
    import glob
    import matplotlib.dates as mdates
    from astropy import units as u
    import matplotlib.pyplot as plt

    lw=2.5
    tickfont=12

    path_to_dodem = '/Users/jmdunca2/do-dem/'



    flare_res = get_saved_flares(add_stdv_flares=True, add_manual_flares=True)
    early_starts = flare_res[0]
    late_stops = flare_res[1]
    
    with open(file, 'rb') as f:
        data = pickle.load(f)
    
    ARDict = data[key]
    
    id_dirs = ARDict['datapaths']
    obsids = ARDict['obsids']
    working_dir = ARDict['working_dir']
    # prepped_aia_dir = ARDict['prepped_aia']
    method = ARDict['method']
    
    # if method=='double':
    #     gauss_stats = ARDict['gauss_stats']
    #     sep_axis = gauss_stats[0][0]
    # else:
    #     sep_axis = ''
        
    
    if method in ['input', 'double']:
        directories = get_region_directories(key, targets_file=file)
        #all_all_time_intervals is a list... 
        #       with entries (lists) for each directory/region
        #             those lists have entries (lists) for each orbit
        #                   those lists have entries for each time interval.
        all_all_time_intervals, fixit = tis.region_time_intervals(directories, id_dirs, shush=True)


        #for each directory
        for d in range(0, len(directories)):
            print(d)
            all_sumcons=[]
            all_sumcons_non=[]
            all_sumcons_flare=[]
            if time_weighted:
                all_sumcons_tw=[]
                all_sumcons_tw_non=[]
                all_sumcons_tw_flare=[]
                
            dir_ = directories[d]
            all_orbits_time_intervals = all_all_time_intervals[d]
            #for each orbit
            for ot in range(0, len(all_orbits_time_intervals)):
                orbittimes = all_orbits_time_intervals[ot] 
                if not orbittimes:
                    continue
                # if key == '26-jul-16_1' and ot==2:
                #     print(orbittimes)
                #     print(orbittimes[-1])
                sumcons = []
                sumcons_non=[]
                sumcons_flare=[]
                if time_weighted:
                    sumcons_tw=[]
                    sumcons_tw_non=[]
                    sumcons_tw_flare=[]
                    
                times = []
                endtime = []
                #for each time interval
                for tt in orbittimes:
                    quitthistime=False
                    first=True
                    timestring = viz.make_timestring(tt)
                    print(timestring)
                    files = glob.glob(dir_+'/'+timestring+'/*'+key+'_MC_'+'*result.pickle')
                    files.sort()
                    #print(files)
                    if not files:
                        print('No DEM Result files found for ', timestring)
                        continue
                    sumcon = 0.5
                    for f in files:
                        #print(f)
                        data, timestring, time = viz.load_DEM(f)
                        if first:
                            first=False
                            res = check_avg_rej(time, data['nustar_datapath'], threshold=accthreshold)
                            if not res[1]:
                                print('For time, ', time[0].strftime('%D %H-%M-%S'), '-', 
                                      time[1].strftime('%D %H-%M-%S'), ' mean accepted events, ', 
                                      res[0], ' below threshold, ', accthreshold)
                                quitthistime=True
                                continue  

                            
                            times.append(time[0].datetime)
                            endtime=time[1]
                            
                        #print(viz.checkresid(data))
                        sumcon += (1 if viz.checkresid(data) == 0 else 0)

                    if quitthistime:
                        continue

                    flare = check_for_flare(time, early_starts, late_stops)
                    
                    if time_weighted:
                        dur = (time[1]-time[0]).to(u.s)
                        time_mult = round(dur.value/seconds_per)
                        sumcons_tw.extend([sumcon for i in range(0,time_mult)])
                        if flare:
                            sumcons_tw_flare.extend([sumcon for i in range(0,time_mult)])
                        else:
                            sumcons_tw_non.extend([sumcon for i in range(0,time_mult)])
                        
                    sumcons.append(sumcon)
                    if flare:
                        sumcons_flare.append(sumcon)
                    else:
                        sumcons_non.append(sumcon)

                if not sumcons:
                    print('No times found in this orbit: ', ot)
                    continue
                    
                if not fetch_for_all:    
                    #times = [t[0].datetime for t in orbittimes]  
                    #times.append(orbittimes[-1][1].datetime)
                    times.append(endtime.datetime)
                    
                    fig, ax = plt.subplots(1, 1, figsize=(15,4), tight_layout = {'pad': 1})
                    ax.stairs(sumcons, times, color='Purple',linewidth=lw, baseline=None)
                    
                    ax.set_yticks([0.5,1.5,2.5,3.5, 4.5, 5.5], ['5', '7', '10', '12.6', '15.8', '(>15.8)'])
                    ax.set_ylabel('Upper Limit (MK)', fontsize=15, color='Black')
                    ax.set_xlabel('Time (UT; HH:MM:SS)', fontsize=15)
                    #ax.set_xlabel('Time (2018 May 29, UTC)', fontsize=15)
                    ax.set_title('Lowest-temperature upper bound where DEM is consistent with NuSTAR 6-10 keV emission',
                                 fontsize=15)
                    ax.grid(True,which='major', axis='y', lw=0.5,color='black', linestyle='dotted')                 
                    #ax.grid(True,which='major', axis='x', lw=0.5,color='black', linestyle='dotted') 
    
                    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                    ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
                    ax.tick_params(axis='x', labelsize=tickfont)
                    ax.tick_params(axis='y', labelsize=tickfont)
                    #ax.legend(loc='upper right', fontsize=15)
            
                    plt.savefig(dir_+'/'+obsids[ot]+'_one_orbit_consistency.png')
                    if not show:
                        plt.close()   

                all_sumcons.extend(sumcons)
                all_sumcons_non.extend(sumcons_non)
                all_sumcons_flare.extend(sumcons_flare)
                if time_weighted:
                    all_sumcons_tw.extend(sumcons_tw)
                    all_sumcons_tw_non.extend(sumcons_tw_non)
                    all_sumcons_tw_flare.extend(sumcons_tw_flare)
                    if fetch_for_all:
                        return all_sumcons_tw, all_sumcons_tw_non, all_sumcons_tw_flare

            if fetch_for_all:
                return all_sumcons, all_sumcons_non, all_sumcons_flare

            if return_region_quiet_max:
                if all_sumcons_tw_non:
                    return np.max(all_sumcons_tw_non)
                else:
                    return
                
            
            if time_weighted:
                if len(all_sumcons_tw) == 0:
                    return
                consistency_hist_plot(all_sumcons_tw, all_sumcons_tw_non, all_sumcons_tw_flare, dir_, label_=key+'_region_'+str(d), 
                               time_weighted=True, ylabel='# of '+str(seconds_per)+'s intervals', show=show)
            else:
                if len(all_sumcons) == 0:
                    return
                consistency_hist_plot(all_sumcons, all_sumcons_non, all_sumcons_flare, dir_, label_=key+'_region_'+str(d), 
                               time_weighted=False, ylabel='# of intervals', show=show)
            
    
    if method=='fit':
        all_time_intervals, all_time_intervals_list = tis.find_all_intervals(working_dir, shush=True)

        all_sumcons=[]
        all_sumcons_non=[]
        all_sumcons_flare=[]
        if time_weighted:
            all_sumcons_tw=[]
            all_sumcons_tw_non=[]
            all_sumcons_tw_flare=[]
            
        #for each orbit
        for ot in range(0, len(all_time_intervals)):
            orbittimes = all_time_intervals[ot] 
            sumcons = []
            sumcons_non=[]
            sumcons_flare=[]
            if time_weighted:
                sumcons_tw=[]
                sumcons_tw_non=[]
                sumcons_tw_flare=[]
            times = []
            endtime = []
            #for each time interval
            for tt in orbittimes:
                first=True
                timestring = viz.make_timestring(tt)
                print(timestring)
                #files = glob.glob(working_dir+'/'+timestring+'/'+'*result.pickle')
                files = glob.glob(working_dir+'/'+timestring+'/*'+key+'_MC_'+'*result.pickle')
                files.sort()
                if not files:
                    print('No DEM Result files found for ', timestring)
                    continue
                sumcon = 0.5
                for f in files:
                    #print(f)
                    data, timestring, time = viz.load_DEM(f)
                    if first:
                        times.append(time[0].datetime)
                        endtime=time[1]
                        first=False
                    #print(viz.checkresid(data))
                    sumcon += (1 if viz.checkresid(data) == 0 else 0)

                flare = check_for_flare(time, early_starts, late_stops)

                if time_weighted:
                    dur = (time[1]-time[0]).to(u.s)
                    time_mult = round(dur.value/seconds_per)
                    sumcons_tw.extend([sumcon for i in range(0,time_mult)])
                    if flare:
                        sumcons_tw_flare.extend([sumcon for i in range(0,time_mult)])
                    else:
                        sumcons_tw_non.extend([sumcon for i in range(0,time_mult)])
                    
                sumcons.append(sumcon)
                if flare:
                        sumcons_flare.append(sumcon)
                else:
                    sumcons_non.append(sumcon)

            if not sumcons:
                print('No times found in this orbit: ', ot)
                continue

            if not fetch_for_all:    
                #times = [t[0].datetime for t in orbittimes]  
                #times.append(orbittimes[-1][1].datetime)
                times.append(endtime.datetime)
                fig, ax = plt.subplots(1, 1, figsize=(15,4), tight_layout = {'pad': 1})
                ax.stairs(sumcons, times, color='Purple',linewidth=lw, baseline=None)
                
                ax.set_yticks([0.5,1.5,2.5,3.5, 4.5, 5.5], ['5', '7', '10', '12.6', '15.8', '(>15.8)'])
                ax.set_ylabel('Upper Limit (MK)', fontsize=15, color='Black')
                ax.set_xlabel('Time (UT; HH:MM:SS)', fontsize=15)
                #ax.set_xlabel('Time (2018 May 29, UTC)', fontsize=15)
                ax.set_title('Lowest-temperature upper bound where DEM is consistent with NuSTAR 6-10 keV emission',
                             fontsize=15)
                ax.grid(True,which='major', axis='y', lw=0.5,color='black', linestyle='dotted')                 
                #ax.grid(True,which='major', axis='x', lw=0.5,color='black', linestyle='dotted') 
    
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
                ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
                ax.tick_params(axis='x', labelsize=tickfont)
                ax.tick_params(axis='y', labelsize=tickfont)
                #ax.legend(loc='upper right', fontsize=15)
        
                plt.savefig(working_dir+'/'+obsids[ot]+'_one_orbit_consistency.png')
                if not show:
                    plt.close()

            all_sumcons.extend(sumcons)
            all_sumcons_non.extend(sumcons_non)
            all_sumcons_flare.extend(sumcons_flare)
            if time_weighted:
                all_sumcons_tw.extend(sumcons_tw)
                all_sumcons_tw_non.extend(sumcons_tw_non)
                all_sumcons_tw_flare.extend(sumcons_tw_flare)
                if fetch_for_all:
                        return all_sumcons_tw, all_sumcons_tw_non, all_sumcons_tw_flare

            if fetch_for_all:
                return all_sumcons, all_sumcons_non, all_sumcons_flare

        if return_region_quiet_max:
            if all_sumcons_tw_non:
                return np.max(all_sumcons_tw_non)
            else:
                return
                

        if time_weighted:
            if len(all_sumcons_tw) == 0:
                return
            consistency_hist_plot(all_sumcons_tw, all_sumcons_tw_non, all_sumcons_tw_flare, working_dir, label_=key, 
                               time_weighted=True, ylabel='# of '+str(seconds_per)+'s intervals', show=show)
        else:
            if len(all_sumcons) == 0:
                return
            consistency_hist_plot(all_sumcons, all_sumcons_non, all_sumcons_flare, working_dir, label_=key+'_region_'+str(d), 
                           time_weighted=False, ylabel='# of intervals', show=show)
            
    

def dem_respects_min_loci(data, plot=True):
    
    #print(data['dn_in'])

    all_loci = np.zeros((len(data['chanax']), len(data['ts_'])))
    for i in np.arange(len(data['chanax'])):
        all_loci[i,:] =  data['dn_in'][i]/data['trmatrix'][:,i]
    
    min_loci = np.min(all_loci, axis=0)
    min_loci_interp = np.interp(data['ts'],data['ts_'],min_loci)

    if plot:
        fig = plt.figure(figsize=(5, 4))
        plt.semilogy(data['ts_'], min_loci)
        plt.semilogy(data['ts'], data['DEM'])
        plt.semilogy(data['ts'], min_loci_interp)
    
        plt.xlim(5.7,7.2)
        plt.ylim(1e20,1e30)


    ahhh = np.where(data['DEM'] > min_loci_interp)[0]

    if len(ahhh) > 0:
        return False
    else:
        return True


def check_loci(key, file, searchstring='_no_xrt_'):

    """
    file: contains targets directory (indexed by key)
    """


    violations = []
    

    
    res_files, tworegions = get_key_resultfiles(key, file, fromhome=True, shush=True,
                        withparams=False,
                        namesearchstring=key+searchstring)    

    
    if tworegions:
        for j in range(0,2):
            #label = ARDict['NOAA_ARID'][j]+'-'+key
            resfiles = res_files[j]
            for r in resfiles:
                #print(r)
                data, timestring, time = viz.load_DEM(r)
                if not dem_respects_min_loci(data, plot=False):
                    print('AHHHHHHHHHHHHHHHH')
                    violations.append(r)
                    



    else:
        for r in res_files:
            data, timestring, time = viz.load_DEM(r)
            if not dem_respects_min_loci(data, plot=False):
                    print('AHHHHHHHHHHHHHHHH')
                    violations.append(r)



    return violations



def check_file_instruments_and_flare(r, early_starts, late_stops, lenrange=[6,6], 
                                     checkacc=True, accthreshold=95, shush=False):

    data, timestring, time = viz.load_DEM(r)
    acc=True
    if checkacc:
        res = check_avg_rej(time, data['nustar_datapath'], threshold=accthreshold)
        if not res[1]:
            if not shush:
                print('For time, ', time[0].strftime('%D %H-%M-%S'), '-', 
                      time[1].strftime('%D %H-%M-%S'), ' mean accepted events, ', 
                      res[0], ' below threshold, ', accthreshold)
    
            acc=False

    res = viz.get_DEM_params(r, save_params_file=True)
    if not res:
        return

    m1, max1, above5_, above7_, above10_, \
        above_peak, below_peak, above_635, below_635, \
           chanax, dn_in, edn_in, \
                powerlaws, EMT_all, EMT_thresh = res


    if np.logical_and(len(chanax) >= lenrange[0], len(chanax) <= lenrange[1]):
        
        data, timestring, time = viz.load_DEM(r)
        flare = check_for_flare(time, early_starts, late_stops)
        #print(flare)
        return above10_, flare, time, acc

    else:
        #print('Length of channels list not in range: ', lenrange)
        #print(chanax)
        return


def sorted_resfiles_dict(file, checkacc=True, accthreshold=95):


    with open(file, 'rb') as f:
        keydict = pickle.load(f)

    keys = keydict.keys()

    dictz = {'all regions': {'flare files': [],
                           'quiet files': []}
            }

    early_starts, late_stops = get_saved_flares(flarepath='./reference_files/', 
                                                add_stdv_flares=True, add_manual_flares=True)

    
    
    for key in keys:
        conditions = ['onlyaia', 'aiaxrt', 'no_xrt', key+'_MC']
        #conditions = ['onlyaia', 'no_xrt']
        lenranges = [[6,6], [7,9], [9,9], [9,12]]
    
        ardict = keydict[key]
        dictz[key+' region_0'] = {'key': key}
        if len(ardict['loc']) == 2:
            dictz[key+' region_1'] = {'key': key}
            
            
        
        for c in range(0, len(conditions)):
            
            cc = conditions[c]               
            res_files_oa, tworegions = get_key_resultfiles(key, file, fromhome=True, 
                                withparams=False,
                                namesearchstring=cc,
                                shush=True)
    
            if c == 3:
                cc = 'all-inst'
            
            if tworegions:
                for j in range(0,2):
                    flarefiles = []
                    quietfiles = []
                    flaretimes = []
                    quiettimes = []
                    quietxrttimes = []
                    flarexrttimes = []
                    rejfiles = []
                    rejtimes = []
                    reglab = 'region_'+str(j)
                    #label = ARDict['NOAA_ARID'][j]+'-'+key
                    resfiles = res_files_oa[j]
                    for r in resfiles:
                        #print(r)
                        res = check_file_instruments_and_flare(r, early_starts, late_stops, 
                                                               checkacc=checkacc, accthreshold=accthreshold, 
                                                               lenrange=lenranges[c])
                        if res is not None:
                            a10, flare, time, acc = res
                            if acc:
                                if flare:
                                    flarefiles.append(r.split('.p')[-2]+'_withparams.pickle')
                                    flaretimes.append(time)
                                else:
                                    #print(r)
                                    quietfiles.append(r.split('.p')[-2]+'_withparams.pickle')
                                    quiettimes.append(time)
        
                                if c == 1:
                                    if flare:
                                        flarexrttimes.append(time)
                                    else:
                                        quietxrttimes.append(time)
                            else:
                                rejfiles.append(r.split('.p')[-2]+'_withparams.pickle')
                                rejtimes.append(time)
                        
        
                    dictz[key+' '+reglab]['flare files '+cc] = flarefiles
                    dictz[key+' '+reglab]['quiet files '+cc] = quietfiles
                    dictz[key+' '+reglab]['flare times '] = flaretimes
                    dictz[key+' '+reglab]['quiet times '] = quiettimes
                    dictz[key+' '+reglab]['rejected times'] = rejtimes
                    dictz[key+' '+reglab]['rejected files'] = rejfiles                 
                    if c == 1:
                        dictz[key+' '+reglab]['flare xrt times '] = flarexrttimes
                        dictz[key+' '+reglab]['quiet xrt times '] = quietxrttimes
    
                    if c == 3:
                        dictz['all regions']['flare files'].extend(flarefiles)
                        dictz['all regions']['quiet files'].extend(quietfiles)
                                             
            
            else:
                flarefiles = []
                quietfiles = []
                flaretimes = []
                quiettimes = []
                quietxrttimes = []
                flarexrttimes = []
                rejfiles = []
                rejtimes = []
                resfiles = res_files_oa
                reglab = 'region_0'
                for r in resfiles:
                    #print(r)
                    res = check_file_instruments_and_flare(r, early_starts, late_stops, 
                                                           checkacc=checkacc, accthreshold=accthreshold, 
                                                           lenrange=lenranges[c])
                    if res is not None:
                        a10, flare, time, acc = res
                        if acc:
                            if flare:
                                flarefiles.append(r.split('.p')[-2]+'_withparams.pickle')
                                flaretimes.append(time)
                            else:
                                #print(r)
                                quietfiles.append(r.split('.p')[-2]+'_withparams.pickle')
                                quiettimes.append(time)
        
                            if c == 1:
                                if flare:
                                    flarexrttimes.append(time)
                                else:
                                    quietxrttimes.append(time)

                        else:
                            rejfiles.append(r.split('.p')[-2]+'_withparams.pickle')
                            rejtimes.append(time)
                    
        
                dictz[key+' '+reglab]['flare files '+cc] = flarefiles
                dictz[key+' '+reglab]['quiet files '+cc] = quietfiles
                dictz[key+' '+reglab]['flare times '] = flaretimes
                dictz[key+' '+reglab]['quiet times '] = quiettimes
                dictz[key+' '+reglab]['rejected times'] = rejtimes
                dictz[key+' '+reglab]['rejected files'] = rejfiles 
                if c == 1:
                    dictz[key+' '+reglab]['flare xrt times '] = flarexrttimes
                    dictz[key+' '+reglab]['quiet xrt times '] = quietxrttimes
    
                if c == 3:
                        dictz['all regions']['flare files'].extend(flarefiles)
                        dictz['all regions']['quiet files'].extend(quietfiles)

    print(dictz.keys())
    return dictz


def print_stats(dictt):

    print(dictt['key'])
    print('Flare-time only-AIA DEMS: ', len(dictt['flare files onlyaia']))
    print('Quiet-time only-AIA DEMS: ', len(dictt['quiet files onlyaia']))
    print('Flare-time AIA+XRT DEMS: ', len(dictt['flare files aiaxrt']))
    print('Quiet-time AIA+XRT DEMS: ', len(dictt['quiet files aiaxrt']))    
    print('Flare-time AIA+NuSTAR DEMS: ', len(dictt['flare files no_xrt']))
    print('Quiet-time AIA+NuSTAR DEMS: ', len(dictt['quiet files no_xrt'])) 
    print('Flare-time All-instrument DEMS: ', len(dictt['flare files all-inst']))
    print('Quiet-time All-instrument DEMS: ', len(dictt['quiet files all-inst']))  
    print('Rejected Times: ', len(dictt['rejected files']))
  
    totxrt = len(dictt['flare files aiaxrt']) + len(dictt['quiet files aiaxrt'])
    totany = len(dictt['flare files all-inst']) + len(dictt['quiet files all-inst'])

    #print('flare xrt times ', len(dictt['flare xrt times ']))
    #print('quiet xrt times ', len(dictt['quiet xrt times ']))
    #print('')


    return (totany > 0), (totxrt > 0)
    
    